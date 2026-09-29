"""Cell shape detection for adherent cells of arbitrary shape."""
import numpy as np
from scipy import ndimage as ndi
from skimage import feature, filters, measure, morphology, segmentation


def projection(stack, skip_first=2):
    """Mean projection over z. The lowest plane(s) are skipped (glass / offset noise)."""
    return stack[skip_first:].mean(axis=0)


def segment_cells(stack, pixel_um, nuclei=None, sigma_um=0.8, min_area_um2=60.0, seed_um=3.0, skip_first=2):
    """Return (label image, smoothed projection). Label 0 = background, 1..N = cells.

    1. Gaussian-smooth the z-mean projection so the dotty membrane signal becomes a filled blob.
    2. Threshold (Li, robust for dim images), fill holes, drop small debris.
    3. Split touching cells with a watershed on the distance map. Seeds are nuclei (DAPI) when a
       nuclear stack is given, otherwise distance-map maxima. No round-cell assumption is made.
    4. Drop cells touching the image border (their contour would be incomplete).
    """
    img = projection(stack, skip_first)
    g = filters.gaussian(img, sigma=sigma_um / pixel_um)
    mask = g > filters.threshold_li(g)
    mask = ndi.binary_fill_holes(morphology.binary_opening(mask, morphology.disk(2)))
    mask = morphology.remove_small_objects(mask, int(min_area_um2 / pixel_um**2))

    dist = filters.gaussian(ndi.distance_transform_edt(mask), sigma=1.5)
    elevation = -dist
    if nuclei is not None:
        n = filters.gaussian(projection(nuclei, skip_first), sigma=seed_um / pixel_um)
        pk = feature.peak_local_max(n, min_distance=int(seed_um / pixel_um), labels=mask.astype(int),
                                    threshold_abs=filters.threshold_li(n))
        elevation = -n  # cell borders follow the valleys between nuclei
        seeds = np.zeros(mask.shape, int)
        seeds[tuple(pk.T)] = np.arange(1, len(pk) + 1)
        cc = measure.label(mask)  # blobs without any nucleus keep their own seed
        for r in measure.regionprops(cc):
            if not seeds[cc == r.label].any():
                seeds[tuple(np.round(r.centroid).astype(int))] = seeds.max() + 1
    else:
        seeds = measure.label(morphology.h_maxima(dist, h=seed_um / pixel_um / 2))
    labels = segmentation.watershed(elevation, seeds, mask=mask)

    border = np.unique(np.concatenate([labels[0], labels[-1], labels[:, 0], labels[:, -1]]))
    labels[np.isin(labels, border)] = 0
    out = np.zeros_like(labels)
    for new, old in enumerate([l for l in np.unique(labels) if l], start=1):
        out[labels == old] = new
    return out, g
