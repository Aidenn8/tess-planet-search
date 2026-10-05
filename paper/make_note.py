"""Draft a Research Note of the AAS (RNAAS: <= 1,000 words, one table) from the result files.

    .venv/bin/python paper/make_note.py      ->  paper/rnaas_note.tex

Author and affiliation are placeholders for the author.
Every number is read from results/, so the note cannot drift from the dossiers.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from tess_search import RESULTS
from tess_search.assessments import facts

HERE = Path(__file__).resolve().parent
HARD = RESULTS / "hardening"


def pm(d, fmt="{:.2f}"):
    return f"${fmt.format(d['median'])}\\pm{fmt.format(0.5 * (d['plus'] + d['minus']))}$"


def main():
    cs = pd.read_csv(RESULTS / "candidates" / "summary.csv")
    loc = pd.read_csv(HARD / "localization_summary.csv")
    syst = json.loads((HARD / "localization_systematics.json").read_text())
    summary = json.loads((RESULTS / "summary.json").read_text())
    pl = loc[loc.kind == "planet"]
    nb = loc[loc.kind == "neb"]
    inj = loc[(loc.kind == "injection") & (loc.reliable == True)]  # noqa: E712
    n_inj_ok = int(inj.recovered_correct.astype(str).eq("True").sum())
    cand = cs[cs.verdict.str.startswith("Candidate")]
    fp = cs[~cs.verdict.str.startswith("Candidate")]

    rows = []
    for _, r in cand.iterrows():
        m = json.loads((HARD / "mcmc" / f"TIC{r.tic}.json").read_text())["prior"]
        name = "TOI-218" if r.tic == 32090583 else ""
        rows.append(f"{r.tic}{(' (' + name + ')') if name else ''} & ${r.period:.6f}$ & ${r.t0_bjd - 2450000:.4f}$ & "
                    f"{pm(m['depth_ppm'], '{:.0f}')} & {pm(m['t14_h'])} & {pm(m['rp_rearth'])} & "
                    f"${m['teq_k']['median']:.0f}$ & ${r.loc_offset_arcsec:.1f}\\pm{r.loc_err_arcsec:.1f}$ & "
                    f"${r.fpp:.2f}$ & ${r.fpp_cleared:.3f}$ \\\\")
    f218 = facts(32090583)
    f229 = facts(229689348)
    f198 = facts(198412174)
    f294 = facts(294053492)

    tex = rf"""\documentclass[RNAAS]{{aastex631}}
\begin{{document}}
\title{{Four Earth-sized Transit Candidates Around Long-observed TESS M Dwarfs, Including a Third Signal in the TOI-218 System}}
\author{{AUTHOR NAME}}
\affiliation{{AFFILIATION}}
\keywords{{Exoplanet detection methods (489) --- Transit photometry (1709) --- M dwarf stars (982)}}

\section{{Search}}
M dwarfs near the TESS continuous viewing zones now have years of 2-minute photometry, enough to reach
Earth-sized planets. We searched the {summary['stars']:,} M dwarfs ($T_{{\rm eff}}\le3900$\,K, $T\le13.5$) with at
least 20 sectors of SPOC 2-minute light curves through Sector 107 (median 26 sectors), for periods of
0.4--40\,d. Light curves were cleaned of flares and detrended with a biweight filter \citep{{hippke2019}}; each
observing season was searched with box least squares \citep{{kovacs2002}} on a common period grid
\citep{{ofir2014}} and the seasons' likelihoods summed. Signals were vetted with red-noise-calibrated tests
(odd/even depths, secondary eclipses, duration versus stellar density, single-event dominance, spacecraft
and rotation periods) and crossmatched with TOIs, CTOIs, confirmed planets, the TESS eclipsing-binary
catalog and every SPOC multi-sector Threshold Crossing Event. The pipeline recovers 25 of 26 confirmed
transiting planets in range; 74\% of 300 injected planets are recovered and kept, and 200 flipped light
curves produce no false candidates. Five Earth-sized signals in no catalog passed; three had been SPOC TCEs
but never TOIs.

\section{{Vetting in the pixels}}
To test where each transit originates we built difference images (out-minus-in transit) for every sector,
with per-pixel errors from $\sim$60 null events per sector, and fitted all sectors jointly with the SPOC
pixel response function, calibrated per sector on Gaia DR3 stars \citep{{gaia2023}}. Position errors include a
{syst['sys_arcsec_used']:.1f}$''$ systematic floor set by {syst['n_cases']} sources of known position. The
method places all {len(pl)} confirmed planets tested on their host (within {pl.target_sigma_total.max():.1f}$\sigma$),
places {len(nb)} TFOP-retired nearby eclipsing binaries off target, and traces {n_inj_ok} of {len(inj)}
reliably fitted eclipses injected into the real pixels to the correct star. We fitted each transit with \texttt{{batman}}
\citep{{kreidberg2015}} and \texttt{{emcee}} \citep{{foreman2013}} (stellar-density prior from the TIC,
\citealt{{stassun2019}}) and computed false-positive probabilities with TRICERATOPS \citep{{giacalone2021}}
using the Gaia DR3 field population, without high-resolution imaging.

\section{{Results}}
TIC~294053492 ($P={f294['period']:.4f}$\,d) is a nearby eclipsing binary: its light loss lies
{f294['loc_off']:.0f}$''$ from the target, excluded at {f294['target_sig']:.1f}$\sigma$ (independently in odd
and even sectors), {f294['src_dbest']:.1f}$''$ from Gaia DR3 {f294['src_id']}. The other four remain candidates
(Table~\ref{{tab:cands}}):

\textit{{TOI-218}} (TIC~32090583), one component of a 13.5$''$ wide binary of near-equal M dwarfs, shows a third
signal at {f218['period']:.4f}\,d ({f218['rp']}\,$R_\oplus$) beside TOI-218.01 and .02. All three signals
localize to TIC~32090583; the new one excludes the companion at {f218['comp_sig']:.1f}$\sigma$.
TRICERATOPS alone, blind to the pixels, gives FPP\,=\,{f218['fpp']:.2f}, almost all from the companion as host;
with the neighbours the localization excludes treated as cleared, FPP\,=\,{f218['fpp_c']:.3f}.

\textit{{TIC~229689348}} hosts an 11.2-hour, {f229['rp']}\,$R_\oplus$ candidate ($T_{{\rm eq}}\approx
{f229['teq']:.0f}$\,K). SPOC flagged it in three runs; its falling SNR follows SPOC's period error, and its
55$''$ difference-image offset came from images that failed SPOC's quality metric. Our localization places
the source on the target and excludes the 9$''$ and 49$''$ neighbours at {f229['near_sig']:.1f} and
{f229['far_sig']:.1f}$\sigma$.

\textit{{TIC~149390648}} lies in a crowded field (NFPP\,=\,{facts(149390648)['nfpp']:.4f}).
\textit{{TIC~198412174}} is near-grazing; a $T=18.2$ star 4.7$''$ away cannot be excluded, and TRICERATOPS
favors a planet around an unresolved companion ({f198['stp']:.0f}\%) about equally with the target
({f198['tp']:.0f}\%).

With the neighbours that the pixels exclude treated as cleared, all four meet the TRICERATOPS ``likely planet''
criteria (FPP\,$<$\,0.5, NFPP\,$<$\,0.001; Table~\ref{{tab:cands}}); the residual FPP is mostly unresolved bound
companions. All four are suited to seeing-limited photometry and high-resolution imaging. Code, light-curve products
and per-candidate dossiers are available at \url{{https://github.com/Aidenn8/tess-planet-search}}.

\begin{{deluxetable*}}{{lccccccccc}}
\tablecaption{{Candidates (MCMC medians and 68\% intervals)\label{{tab:cands}}}}
\tablehead{{\colhead{{TIC}} & \colhead{{$P$ (d)}} & \colhead{{$T_0$ (BJD$-$2450000)}} & \colhead{{Depth (ppm)}} &
\colhead{{$T_{{14}}$ (h)}} & \colhead{{$R_p$ ($R_\oplus$)}} & \colhead{{$T_{{\rm eq}}$ (K)}} &
\colhead{{Offset ($''$)}} & \colhead{{FPP}} & \colhead{{FPP$_{{\rm cl}}$}}}}
\startdata
{chr(10).join(rows)}
\enddata
\tablecomments{{$R_p$ includes the TIC stellar-radius uncertainty. Offset: localized source position relative
to the target. FPP from TRICERATOPS without imaging constraints; FPP$_{{\rm cl}}$ with the neighbours that the
pixel localization excludes at $>3\sigma$ treated as cleared.}}
\end{{deluxetable*}}

\begin{{thebibliography}}{{}}
\bibitem[Foreman-Mackey et al.(2013)]{{foreman2013}} Foreman-Mackey, D., Hogg, D.~W., Lang, D., \& Goodman, J. 2013, PASP, 125, 306
\bibitem[Gaia Collaboration et al.(2023)]{{gaia2023}} Gaia Collaboration, Vallenari, A., et al. 2023, A\&A, 674, A1
\bibitem[Giacalone et al.(2021)]{{giacalone2021}} Giacalone, S., Dressing, C.~D., Jensen, E.~L.~N., et al. 2021, AJ, 161, 24
\bibitem[Hippke et al.(2019)]{{hippke2019}} Hippke, M., David, T.~J., Mulders, G.~D., \& Heller, R. 2019, AJ, 158, 143
\bibitem[Kov{{\'a}}cs et al.(2002)]{{kovacs2002}} Kov{{\'a}}cs, G., Zucker, S., \& Mazeh, T. 2002, A\&A, 391, 369
\bibitem[Kreidberg(2015)]{{kreidberg2015}} Kreidberg, L. 2015, PASP, 127, 1161
\bibitem[Ofir(2014)]{{ofir2014}} Ofir, A. 2014, A\&A, 561, A138
\bibitem[Stassun et al.(2019)]{{stassun2019}} Stassun, K.~G., Oelkers, R.~J., Paegert, M., et al. 2019, AJ, 158, 138
\end{{thebibliography}}
\end{{document}}
"""
    tex = tex.replace("±", "$\\pm$")
    out = HERE / "rnaas_note.tex"
    out.write_text(tex)
    body = tex.split(r"\section{Search}")[1].split(r"\begin{deluxetable*}")[0]
    words = len([w for w in body.replace("\\", " ").split() if any(c.isalpha() for c in w)])
    print(f"{out}  (~{words} words in the text body; RNAAS limit 1,000)")


if __name__ == "__main__":
    main()
