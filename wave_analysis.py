#!/usr/bin/env python3
"""
Are Bridget Riley's waves actually sine waves?

Fourier analysis of the wave contours in three paintings by Bridget Riley:
'Study for Arrest Series' (1965), 'Cataract 3' (1967) and 'Gala' (1974).

Pipeline:
  1. Load image, grayscale, Otsu binarize (pigment vs near-white ground).
  2. Connected components -> each painted stripe.
  3. For each stripe, trace its top and bottom boundary as y(x), one sample
     per pixel column.  Each boundary is one "copy" of the wave.
  4. Per copy: estimate the fundamental wavelength (FFT guess +
     least-squares refinement), then fit a Fourier series
        y(x) ~ a0 + a1*x_lin + sum_n cn * cos(n*w*x - phi_n)
     by linear least squares at the refined wavelength.
  5. Phase offsets between copies are removed by reporting amplitudes cn/c1
     (translation-invariant) and aligned phases psi_n = phi_n - n*phi_1
     (also translation-invariant).
  6. Aggregate across copies: mean/std spectrum, circular-mean phases,
     THD and RMS-deviation-from-pure-sine metrics.  Plots + printed tables.
  7. Similarity analysis: cross-correlate the copies' deviation-from-sine
     profiles within each painting (template test) and the mean profiles
     between paintings (shared-template test).

Outputs: out/<tag>_analysis.png per painting, out/similarity.png, and
tables on stdout.  See readme.md for the findings.
"""

import numpy as np
import cv2
import matplotlib

matplotlib.use("Agg")  # headless: we only save PNGs, never open a window
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

# Number of harmonics in every fit.  12 is far past where the paintings'
# coherent content dies out (~harmonic 8), so the top harmonics measure the
# noise floor rather than the curve.
NHARM = 12


# ----------------------------------------------------------------------
# contour extraction
# ----------------------------------------------------------------------

def extract_edges(path, min_width_frac=0.90):
    """Return (binary image, list of (x, y) boundary curves).

    y is in image coordinates (down = positive); flipped to 'up = positive'
    before analysis.  Edges clipped by the image frame are discarded.
    Otsu on grayscale separates pigment from ground for both the black/white
    works and the coloured Cataract stripes (any colour vs near-white).
    """
    gray = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    if gray is None:
        raise FileNotFoundError(path)
    H, W = gray.shape
    # THRESH_BINARY_INV: pigment is darker than ground, so invert to make
    # the stripes the foreground (255).  Otsu picks the threshold from the
    # bimodal histogram automatically.
    _, bw = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    ink = (bw > 0)  # True where pigment

    # Each stripe becomes one connected component.
    nlab, labels = cv2.connectedComponents(ink.astype(np.uint8))
    edges = []
    for lab in range(1, nlab):  # label 0 is the background
        mask = labels == lab
        cols = mask.any(axis=0)
        xs = np.where(cols)[0]
        if len(xs) < min_width_frac * W:
            continue  # speck or partial blob, not a full stripe
        x0, x1 = xs[0], xs[-1]
        sub = mask[:, x0:x1 + 1]
        # Per column: first foreground pixel from the top and from the
        # bottom.  argmax on a boolean column returns the first True.
        top = sub.argmax(axis=0)
        bot = sub.shape[0] - 1 - sub[::-1].argmax(axis=0)
        colsum = sub.sum(axis=0)
        # Keep only columns where the stripe is a single solid run; columns
        # with holes or merged neighbours would give a false top/bottom.
        solid = (bot - top + 1) == colsum
        x = np.arange(x0, x1 + 1)
        # The top and bottom boundary of a stripe are two independent
        # copies of the wave.
        for y in (top.astype(float), bot.astype(float)):
            good = solid & (y > 0) & (y < H - 1)  # drop frame-clipped pixels
            if good.mean() < 0.97:
                continue  # too damaged to be a trustworthy waveform
            edges.append((x[good].astype(float), y[good]))
    return ink, edges


# ----------------------------------------------------------------------
# per-copy Fourier fit
# ----------------------------------------------------------------------

def harmonic_design(x, lam, nharm):
    """Least-squares design matrix: [1, x, cos(nwx), sin(nwx) for n=1..N].

    The DC and linear columns absorb the stripe's vertical offset and any
    slow drift/tilt, so they do not leak into the harmonic estimates.
    """
    w = 2 * np.pi / lam
    cols = [np.ones_like(x), x - x.mean()]  # DC + linear drift
    for n in range(1, nharm + 1):
        cols.append(np.cos(n * w * x))
        cols.append(np.sin(n * w * x))
    return np.column_stack(cols)


def fit_at_lambda(x, y, lam, nharm):
    """Linear least-squares Fourier fit at a fixed wavelength.

    Fitting by least squares instead of taking an FFT handles a non-integer
    number of periods without windowing/leakage artefacts.  Returns
    (coefficients, RMS residual).
    """
    A = harmonic_design(x, lam, nharm)
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    return coef, float(np.sqrt(np.mean(resid ** 2)))


def estimate_lambda(x, y):
    """FFT-based initial wavelength guess, refined by least squares.

    The FFT (Hann window, 8x zero-padding for peak interpolation) gets
    within a few percent; the scalar minimizer then tunes lambda to
    minimize the residual of a 3-harmonic fit, which is insensitive to
    spectral leakage.
    """
    # Detrend before the FFT so the linear ramp doesn't swamp the spectrum.
    y0 = y - np.polyval(np.polyfit(x, y, 1), x)
    n = len(y0)
    pad = 1 << (int(np.ceil(np.log2(n))) + 3)  # 8x zero-padding
    spec = np.abs(np.fft.rfft(y0 * np.hanning(n), pad))
    freqs = np.fft.rfftfreq(pad, d=1.0)
    # Ignore frequencies below 1.5 cycles over the span: residual detrend
    # power lives there and would steal the peak.
    lo = np.searchsorted(freqs, 1.5 / n)
    k = lo + spec[lo:].argmax()
    lam0 = 1.0 / freqs[k]
    # Refine within +/-15% of the FFT guess; 3 harmonics are enough to lock
    # onto the fundamental without chasing noise.
    r = minimize_scalar(lambda lam: fit_at_lambda(x, y, lam, 3)[1],
                        bounds=(0.85 * lam0, 1.15 * lam0), method="bounded")
    return float(r.x)


def analyse_edge(x, y_img, nharm=NHARM):
    """Fit one boundary curve.  Returns dict or None if unusable."""
    y = -y_img  # image y points down; flip so 'up' is positive
    lam = estimate_lambda(x, y)
    nper = (x[-1] - x[0]) / lam
    if nper < 2.0:
        return None  # under two periods: wavelength/harmonics ill-determined
    coef, rms = fit_at_lambda(x, y, lam, nharm)
    # coef layout: [DC, linear, a1, b1, a2, b2, ...]
    a = coef[2::2]
    b = coef[3::2]
    # Convert (a, b) quadrature pairs to amplitude/phase form:
    # a*cos + b*sin = c*cos(nwx - phi) with c = |(a,b)|, phi = atan2(b, a).
    c = np.hypot(a, b)
    phi = np.arctan2(b, a)              # y = sum c_n cos(n w x - phi_n)
    # Translation-invariant phases: shifting x by D adds n*w*D to phi_n, so
    # psi_n = phi_n - n*phi_1 cancels the shift.  psi is what two copies of
    # the same curve share regardless of where they sit in the painting.
    psi = (phi - np.arange(1, nharm + 1) * phi[0]) % (2 * np.pi)
    # Pure-sine comparison: same fit with the fundamental only.  Its
    # residual is the total deviation from the best possible sine.
    _, rms_sine = fit_at_lambda(x, y, lam, 1)
    return dict(x=x, y=y, lam=lam, nper=nper, c=c, psi=psi, phi1=phi[0],
                rms_full=rms, rms_sine=rms_sine,
                # total harmonic distortion: harmonic power above the
                # fundamental, as a fraction of the fundamental
                thd=np.sqrt((c[1:] ** 2).sum()) / c[0])


# ----------------------------------------------------------------------
# aggregation and reporting
# ----------------------------------------------------------------------

def circ_mean(ang, w):
    """Weighted circular mean of angles.

    Returns (mean angle, concentration).  Concentration is the resultant
    length in [0, 1]: 1 = all copies agree on the phase (a systematic
    feature of the curve), 0 = phases are random (noise).
    """
    z = np.sum(w * np.exp(1j * ang))
    return np.angle(z), abs(z) / np.sum(w)


def analyse_image(path, tag, outdir="out", min_width_frac=0.90):
    """Full analysis of one painting: fit every copy, print the canonical
    Fourier series and sine-deviation metrics, save the figure.  Returns
    the per-copy fit results (consumed by similarity_analysis)."""
    ink, edges = extract_edges(path, min_width_frac)
    results = [r for r in (analyse_edge(x, y) for x, y in edges) if r]
    if not results:
        raise RuntimeError(f"no usable wave edges found in {path}")

    nharm = NHARM
    C = np.array([r["c"] for r in results])          # (copies, nharm)
    rel = C / C[:, :1]                                # cn / c1 per copy
    rel_mean, rel_std = rel.mean(axis=0), rel.std(axis=0)
    lam_mean = np.mean([r["lam"] for r in results])
    c1_mean = C[:, 0].mean()
    thd = np.array([r["thd"] for r in results])
    rms_sine = np.array([r["rms_sine"] for r in results])
    rms_full = np.array([r["rms_full"] for r in results])

    # Aligned phases psi_n averaged across copies.  Circular mean (angles
    # wrap) weighted by each copy's harmonic amplitude, so copies in which
    # the harmonic is well measured count more.  The concentration column
    # is the key significance test: amplitudes alone are biased upward by
    # noise, but phase agreement across copies cannot be faked by noise.
    psi_mean = np.zeros(nharm)
    psi_conc = np.zeros(nharm)
    for n in range(nharm):
        psi_mean[n], psi_conc[n] = circ_mean(
            np.array([r["psi"][n] for r in results]), C[:, n])
    psi_mean[0] = 0.0  # by construction (psi_1 = phi_1 - 1*phi_1 = 0)

    # ---------------- printed report ----------------
    print(f"\n{'=' * 74}\n{tag}: {path}")
    print(f"  usable wave copies (stripe edges): {len(results)}")
    print(f"  mean wavelength: {lam_mean:7.1f} px   "
          f"mean fundamental amplitude c1: {c1_mean:6.1f} px")
    print(f"  periods per copy: "
          f"{min(r['nper'] for r in results):.1f}-{max(r['nper'] for r in results):.1f}")
    print(f"\n  Fourier series (canonical, fundamental phase = 0):")
    print(f"     y(t) = sum_n  c_n * cos(n*t - psi_n),   t = 2*pi*x/lambda")
    print(f"  {'n':>3} {'c_n/c_1':>9} {'+/-':>7} {'c_n/c_1 dB':>11} "
          f"{'psi_n deg':>10} {'phase coherence':>16}")
    for n in range(nharm):
        db = 20 * np.log10(max(rel_mean[n], 1e-12))
        print(f"  {n + 1:>3} {rel_mean[n]:>9.4f} {rel_std[n]:>7.4f} "
              f"{db:>10.1f}  {np.degrees(psi_mean[n]):>9.1f} {psi_conc[n]:>15.2f}")
    print(f"\n  How close to a pure sine?")
    print(f"    THD (harmonics 2-{nharm}):            "
          f"{100 * thd.mean():.2f}%  (per copy {100 * thd.min():.2f}-{100 * thd.max():.2f}%)")
    print(f"    RMS deviation from best sine:     "
          f"{rms_sine.mean():.2f} px = {100 * rms_sine.mean() / c1_mean:.2f}% of amplitude")
    print(f"    RMS residual of {nharm}-harmonic fit: "
          f"{rms_full.mean():.2f} px = {100 * rms_full.mean() / c1_mean:.2f}% of amplitude"
          f"  (noise floor)")

    # ---------------- figure ----------------
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.4, 1, 1])

    # Panel 1 (full width): binarized painting with every traced edge
    # overlaid — the visual check that the tracing followed the stripes.
    ax = fig.add_subplot(gs[0, :])
    ax.imshow(ink, cmap="gray_r", interpolation="nearest")
    for x, y in edges:
        ax.plot(x[::4], y[::4], lw=0.6)  # decimate 4x, plotting only
    ax.set_title(f"{tag}: binarized image with traced stripe edges "
                 f"({len(results)} usable wave copies)")
    ax.set_axis_off()

    # Panel 2: all copies overlaid on one period after phase alignment and
    # amplitude normalization, against a pure sine.  Drawn from the fitted
    # series (not raw pixels) so copies with different wavelengths overlay.
    ax = fig.add_subplot(gs[1, 0])
    t = np.linspace(0, 2 * np.pi, 400)
    for r in results:
        # shift so the fundamental is cos(t); evaluate that copy's series
        yy = sum(r["c"][n] * np.cos((n + 1) * t - r["psi"][n])
                 for n in range(nharm)) / r["c"][0]
        ax.plot(t, yy, color="0.6", lw=0.7)
    y_mean = sum(rel_mean[n] * np.cos((n + 1) * t - psi_mean[n])
                 for n in range(nharm))
    ax.plot(t, y_mean, "k", lw=2, label="mean wave (12-harmonic series)")
    ax.plot(t, np.cos(t), "r--", lw=1.5, label="pure sine")
    ax.set_xlabel("phase t (rad)")
    ax.set_ylabel("y / c1")
    ax.set_title("all copies, phase-aligned, vs pure sine")
    ax.legend(fontsize=8)

    # Panel 3: the difference between the mean wave and a pure sine —
    # the painting's systematic "shape signature".
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(t, y_mean - np.cos(t), "k", lw=1.5)
    ax.axhline(0, color="r", ls="--", lw=0.8)
    ax.set_xlabel("phase t (rad)")
    ax.set_ylabel("(mean wave - sine) / c1")
    ax.set_title("deviation of mean wave from pure sine")

    nn = np.arange(1, nharm + 1)
    # Per-harmonic amplitude noise: a single sinusoid estimated by least
    # squares in white noise of std `rms` over N samples has amplitude
    # error ~ rms*sqrt(2/N).  This is the honest comparison line for the
    # spectrum bars (the broadband residual itself would overstate it).
    N_mean = np.mean([len(r["x"]) for r in results])
    sigma_c = rms_full.mean() * np.sqrt(2.0 / N_mean) / c1_mean
    # Colour-code significance: green = the copies agree on this
    # harmonic's phase, so it is a real feature of the curve.
    colors = ["seagreen" if psi_conc[n] > 0.5 else "0.65" for n in range(nharm)]

    # Panels 4 + 5: the amplitude spectrum, linear and in dB.
    ax = fig.add_subplot(gs[2, 0])
    ax.bar(nn, rel_mean, yerr=rel_std, color=colors, capsize=3)
    ax.set_xlabel("harmonic n")
    ax.set_ylabel("c_n / c_1")
    ax.set_title("amplitude spectrum (linear) — green = phase-coherent across copies")
    ax.set_xticks(nn)

    ax = fig.add_subplot(gs[2, 1])
    db = 20 * np.log10(np.maximum(rel_mean, 1e-12))
    ax.bar(nn, db, color=colors)
    floor_db = 20 * np.log10(sigma_c)
    ax.axhline(floor_db, color="r", ls="--", lw=1,
               label=f"per-harmonic noise level ({floor_db:.0f} dB)")
    ax.set_xlabel("harmonic n")
    ax.set_ylabel("dB rel. fundamental")
    ax.set_ylim(min(db.min(), floor_db) - 5, 5)
    ax.set_title("amplitude spectrum (dB)")
    ax.set_xticks(nn)
    ax.legend(fontsize=8)

    fig.suptitle(f"{tag} — wave-contour Fourier analysis", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = f"{outdir}/{tag}_analysis.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  figure: {out}")
    return results


# ----------------------------------------------------------------------
# within-painting curve similarity
# ----------------------------------------------------------------------

def deviation_profiles(results, nharm=NHARM, npts=256):
    """Per copy: its systematic deviation-from-sine over one period,
    reconstructed from harmonics 2..nharm, phase-aligned (fundamental at
    phase 0) and normalized by the fundamental amplitude.  Dimensionless,
    so directly comparable between copies and between paintings.

    The fundamental is excluded on purpose: every copy is ~a sine, so
    including it would make all profiles correlate near 1 and hide the
    differences.  What remains is each copy's "error fingerprint".
    """
    t = np.linspace(0, 2 * np.pi, npts, endpoint=False)
    profs = []
    for r in results:
        rel = r["c"] / r["c"][0]
        d = sum(rel[n] * np.cos((n + 1) * t - r["psi"][n])
                for n in range(1, nharm))  # n=1 -> harmonic 2: skip fundamental
        profs.append(d)
    return np.array(profs)


def similarity_analysis(all_results, outdir="out"):
    """N x N cross-correlation of deviation profiles within each painting.

    Pearson r per pair of copies: +1 = both copies deviate from a pure sine
    in exactly the same way (shared template), 0 = unrelated deviations
    (independent hand wobble).  The distribution of the N(N-1)/2 pairwise
    values is the tightness measure; Pearson r is scale-invariant, so the
    distributions are comparable between paintings.
    """
    tags = list(all_results)
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, len(tags), height_ratios=[1.15, 1])

    print(f"\n{'=' * 74}\nwithin-painting curve similarity "
          f"(pairwise correlation of deviation-from-sine profiles)")
    print(f"  {'painting':>10} {'copies':>7} {'pairs':>7} {'median r':>9} "
          f"{'IQR':>13} {'r > 0.5':>8}")

    dists = {}
    for k, tag in enumerate(tags):
        # Sort copies by vertical position in the painting (top first), so
        # structure in the matrix maps onto position on the canvas.  (In
        # Gala this exposes diagonal ridges every ~41 copies, where the
        # stripe-to-stripe phase shift completes a full 360-degree cycle
        # and canvas-fixed distortions re-register; see readme.)
        results = sorted(all_results[tag], key=lambda r: -r["y"].mean())
        profs = deviation_profiles(results)
        corr = np.corrcoef(profs)
        n = len(results)
        off = corr[np.triu_indices(n, 1)]  # the N(N-1)/2 distinct pairs
        dists[tag] = off
        q1, q2, q3 = np.percentile(off, [25, 50, 75])
        print(f"  {tag:>10} {n:>7} {len(off):>7} {q2:>9.2f} "
              f"  [{q1:>5.2f},{q3:>5.2f}] {100 * (off > 0.5).mean():>7.0f}%")

        ax = fig.add_subplot(gs[0, k])
        im = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r",
                       interpolation="nearest")
        ax.set_title(f"{tag}: {n}x{n} correlation matrix\n"
                     f"(copies ordered top to bottom of painting)",
                     fontsize=10)
        ax.set_xlabel("copy")
        ax.set_ylabel("copy")
        fig.colorbar(im, ax=ax, fraction=0.046)

    # Bottom panel: the headline display — distribution of pairwise r per
    # painting.  Tight against +1 = template; broad / centred on 0 = each
    # copy deviates its own way.
    ax = fig.add_subplot(gs[1, :])
    bins = np.linspace(-1, 1, 81)
    colors = dict(zip(tags, ["tab:orange", "tab:green", "tab:blue"]))
    for tag in tags:
        ax.hist(dists[tag], bins=bins, density=True, histtype="stepfilled",
                alpha=0.35, color=colors[tag], label=f"{tag} "
                f"(median r = {np.median(dists[tag]):.2f})")
        ax.hist(dists[tag], bins=bins, density=True, histtype="step",
                color=colors[tag], lw=1.5)
        ax.axvline(np.median(dists[tag]), color=colors[tag], ls="--", lw=1)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("pairwise correlation r of deviation-from-sine profiles")
    ax.set_ylabel("density")
    ax.set_title("distribution of pairwise correlations: right-packed = "
                 "copies share one template; centred on 0 = independent deviations")
    ax.legend()

    fig.suptitle("How similar are the wave copies within each painting?",
                 fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = f"{outdir}/similarity.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  figure: {out}")

    # Did different paintings share a template?  Correlate the paintings'
    # mean deviation profiles.  Phase alignment leaves exactly one
    # registration freedom between two different paintings: a mirror flip
    # (t -> -t maps a phase-zero cosine onto itself), so report both.
    def pearson(a, b):
        a = a - a.mean()
        b = b - b.mean()
        return float(a @ b / np.sqrt((a @ a) * (b @ b)))

    print(f"\n  cross-painting correlation of mean deviation profiles:")
    means = {t: deviation_profiles(all_results[t]).mean(axis=0) for t in tags}
    for i, ta in enumerate(tags):
        for tb in tags[i + 1:]:
            print(f"  {ta:>10} vs {tb:<10}: r = {pearson(means[ta], means[tb]):+.2f}, "
                  f"mirrored r = {pearson(means[ta], means[tb][::-1]):+.2f}")


if __name__ == "__main__":
    import os
    os.makedirs("out", exist_ok=True)
    all_results = {
        "arrest": analyse_image("images/riley_arrest.jpg", "arrest"),
        "cataract3": analyse_image("images/riley_cataract3.jpg", "cataract3"),
        "gala": analyse_image("images/riley_gala.jpg", "gala"),
    }
    similarity_analysis(all_results)
