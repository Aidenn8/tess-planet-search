"""Read SPOC Data Validation (DV) XML reports: NASA's own vetting of a TCE.

The fields used here are the ones a TFOP vetter looks at first:
  * multi-sector difference-image centroid offset from the TIC position
    (msTicCentroidOffsets: where the dimming source sits, with its uncertainty),
    and how many per-sector difference images passed SPOC's own quality metric
  * odd/even depth comparison, weak secondary eclipse
  * ghost diagnostic (core vs halo aperture correlation: a high halo value means
    the signal is scattered light from a bright star elsewhere on the detector)
  * bootstrap false-alarm probability, whether the transit fit converged
"""
import xml.etree.ElementTree as ET


def _tag(el):
    return el.tag.split("}", 1)[-1]


def _find(el, name):
    for e in el.iter():
        if _tag(e) == name:
            return e
    return None


def _f(el, key="value"):
    try:
        return float(el.attrib[key]) if el is not None else None
    except (KeyError, ValueError):
        return None


def read_dv(path):
    """List of dicts, one per TCE (planetResults) in a DV XML file."""
    root = ET.parse(path).getroot()
    out = []
    for pr in (e for e in root.iter() if _tag(e) == "planetResults"):
        fit = _find(pr, "allTransitsFit")
        params = {}
        if fit is not None:
            mp = _find(fit, "modelParameters")
            for p in (mp if mp is not None else []):
                params[p.attrib.get("name")] = (_f(p), _f(p, "uncertainty"))
        cand = _find(pr, "planetCandidate")
        mq = _find(pr, "msTicCentroidOffsets")
        mqc = _find(pr, "msControlCentroidOffsets")
        qual = _find(pr, "summaryQualityMetric")
        row = {
            "planet": int(pr.attrib.get("planetNumber", 0)),
            "fit_converged": fit is not None and fit.attrib.get("fullConvergence") == "true",
            "fit_snr": _f(fit, "modelFitSnr"),
            "period": params.get("orbitalPeriodDays", (None, None))[0],
            "period_err": params.get("orbitalPeriodDays", (None, None))[1],
            "depth_ppm": params.get("transitDepthPpm", (None, None))[0],
            "planet_radius_re": params.get("planetRadiusEarthRadii", (None, None))[0],
            "mes": _f(cand, "maxMultipleEventSigma") if cand is not None else None,
            "n_transits": _f(cand, "observedTransitCount") if cand is not None else None,
            "oddeven_sig": _f(_find(pr, "oddEvenTransitDepthComparisonStatistic"), "significance"),
            "ghost_core": _f(_find(pr, "coreApertureCorrelationStatistic")),
            "ghost_halo": _f(_find(pr, "haloApertureCorrelationStatistic")),
            "bootstrap_pfa": _f(_find(pr, "bootstrapResults"), "significance"),
        }
        ws = _find(pr, "weakSecondary")
        if ws is not None:
            row["weak_secondary_mes"] = _f(ws, "maxMes")
            row["weak_secondary_phase_days"] = _f(ws, "maxMesPhaseInDays")
        for name, el in (("tic", mq), ("control", mqc)):
            if el is not None:
                off = _find(el, "meanSkyOffset")
                row[f"ms_{name}_offset_arcsec"] = _f(off)
                row[f"ms_{name}_offset_err_arcsec"] = _f(off, "uncertainty")
        if qual is not None:
            row["diff_images_good"] = int(qual.attrib.get("numberOfGoodMetrics", 0))
            row["diff_images_attempted"] = int(qual.attrib.get("numberOfAttempts", 0))
        # per-sector difference images: how many had a usable fit
        dis = [e for e in pr if _tag(e) == "differenceImageResults"]
        row["n_diff_images"] = len(dis)
        out.append(row)
    return out
