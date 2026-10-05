"""Gaia neighbour bookkeeping: 19-digit source IDs must survive exactly."""
import numpy as np
import pandas as pd

from tess_search import followup


def test_neighbour_ids_are_exact(monkeypatch):
    # two real Gaia DR3 IDs that are 1,536 apart but round to nearby floats
    ids = np.array([4667466476687967232, 4667466476687968768, 4667466721501827072], dtype=np.int64)
    fake = pd.DataFrame({"source_id": ids, "ra": [58.42, 58.43, 58.44], "dec": [-68.74] * 3,
                         "phot_g_mean_mag": [13.0, 15.0, 15.5], "phot_rp_mean_mag": [12.0, 14.0, 14.5],
                         "ruwe": [1.1, 1.0, 1.0], "non_single_star": [0, 0, 0], "parallax": [20.0, 1.0, 1.0],
                         "sep_arcsec": [0.1, 30.0, 50.0]})
    monkeypatch.setattr(followup, "gaia_neighbours", lambda ra, dec: fake)
    out = followup.neighbour_check(58.42, -68.74, 1e-3, target_gaia_id=str(ids[0]))
    got = [n["source_id"] for n in out["neighbours"]]
    assert got == ids[1:].tolist()
    assert len(set(got)) == len(got)
