"""Receptor clustering metric: coefficient of variation (CV) of membrane intensity.

CV = SD / (mean - background) of the intensity sampled along the cell outline.
Evenly spread receptor -> flat profile -> low CV; receptor patches/aggregates -> peaky profile -> high CV.
Dividing by the mean makes CV independent of expression level.
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import filters, measure


def estimate_background(plane, labels, pixel_um, margin_um=2.0):
    """Median of all pixels farther than `margin_um` from any detected cell."""
    far = ~ndi.binary_dilation(labels > 0, iterations=max(1, int(margin_um / pixel_um)))
    return float(np.median(plane[far]))


def best_plane(stack, cell_mask, pixel_um, reference=None, band_um=1.0, skip_first=2):
    """z-index of the cell equator.

    With a nuclear stack the equator is where the nucleus is brightest (adherent cells: the bottom planes
    only show the glass-contact footprint and the top planes the cap). Without it, the plane with the
    brightest outline band is used. The first `skip_first` planes (glass glare) are never picked.
    """
    if reference is not None:
        score = reference[:, cell_mask].mean(axis=1)
    else:
        n = max(1, int(band_um / pixel_um))
        ring = ndi.binary_dilation(cell_mask, iterations=n) & ~ndi.binary_erosion(cell_mask, iterations=n)
        score = stack[:, ring].mean(axis=1)
    score = ndi.uniform_filter1d(score, 3)
    score[:skip_first] = -np.inf
    return int(np.argmax(score))


def equatorial_plane(stack, z, half_width=1):
    """Mean of the planes z-half_width..z+half_width (less shot noise; ~1 depth of field thick)."""
    return stack[max(0, z - half_width): z + half_width + 1].mean(axis=0)


def membrane_profile(plane, cell_mask, pixel_um, band_um=0.8, smooth_px=1.0):
    """Sample the intensity along the cell outline.

    The outline is traced on the (arbitrary-shaped) cell mask, resampled at 1 px spacing, and at every
    outline point the intensity is averaged across a band of +-band_um perpendicular to the outline.
    Returns (distance_um[N], intensity[N], outline_yx[N, 2]).
    """
    smooth_mask = filters.gaussian(cell_mask.astype(float), sigma=1.5)
    contour = max(measure.find_contours(smooth_mask, 0.5), key=len)
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(contour, axis=0).T))]
    t = np.arange(0, seg[-1], 1.0)
    pts = np.c_[np.interp(t, seg, contour[:, 0]), np.interp(t, seg, contour[:, 1])]

    tangent = np.gradient(ndi.uniform_filter1d(pts, 5, axis=0, mode="wrap"), axis=0)
    normal = np.c_[tangent[:, 1], -tangent[:, 0]]
    normal /= np.linalg.norm(normal, axis=1, keepdims=True) + 1e-9

    img = filters.gaussian(plane, sigma=smooth_px, preserve_range=True)
    offsets = np.arange(-band_um / pixel_um, band_um / pixel_um + 1e-6, 0.5)
    samples = [ndi.map_coordinates(img, (pts + o * normal).T, order=1, mode="nearest") for o in offsets]
    return t * pixel_um, np.mean(samples, axis=0), pts


def cell_metrics(profile, background):
    mean, sd = float(np.mean(profile)), float(np.std(profile, ddof=1))
    signal = mean - background
    return dict(Membrane_Mean=mean, SD=sd, Background=background, Signal=signal,
                CV=sd / signal if signal > 0 else np.nan)
