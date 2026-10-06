"""Build the project website (GitHub Pages) from the repository's results and write-up.

    .venv/bin/python site/build.py          ->  site/_site/

Pages:
    index.html             research-style overview (numbers read from results/)
    writeup/index.html     chapter list
    writeup/NN-*.html      docs/NN-*.md rendered with MathJax
    follow-up.html         FOLLOW_UP.md plus upcoming observable transits
The output directory is published on the gh-pages branch.
"""
import csv
import html
import json
import re
import shutil
import sys
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "site" / "_site"
REPO = "https://github.com/Aidenn8/tess-planet-search"
BRANCH = "research"
SITE = "https://aidenn8.github.io/tess-planet-search"
AUTHOR = "Aidenn8"                      # display name on the site and in the citation
AUTHOR_URL = "https://github.com/Aidenn8"
TITLE = "Four Earth-sized planet candidates from a deep search of the longest-observed TESS M dwarfs"
DESCRIPTION = ("A transit search of the 1,279 M dwarfs NASA's TESS has observed the longest, with every "
               "candidate checked in the light curve, in the pixels, statistically, and against NASA's own pipeline.")

sys.path.insert(0, str(ROOT))

ICON = {
    "code": '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="m16 18 6-6-6-6M8 6l-6 6 6 6"/></svg>',
    "doc": '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/><path d="M14 2v6h6M8 13h8M8 17h6"/></svg>',
    "list": '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M8 6h13M8 12h13M8 18h13M3 6h.01M3 12h.01M3 18h.01"/></svg>',
    "scope": '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="9"/><circle cx="12" cy="12" r="3"/><path d="M12 3v3M12 18v3M3 12h3M18 12h3"/></svg>',
    "table": '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="3" y="4" width="18" height="16" rx="2"/><path d="M3 10h18M3 15h18M9 4v16"/></svg>',
}
FAVICON = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32"><circle cx="16" cy="16" r="13" fill="#eb6834"/>'
           '<circle cx="11" cy="13" r="4" fill="#111110"/></svg>')

CHAPTERS = [
    ("01-data-and-search", "Data and search", "sample, cleaning, stacked seasonal BLS, crossmatch"),
    ("02-vetting", "Vetting", "red-noise-aware tests, thresholds, verdicts"),
    ("03-sensitivity-and-reliability", "Sensitivity and reliability", "known planets, 300 injections, 200 inverted light curves"),
    ("04-pixel-level-localization", "Pixel-level localization", "difference imaging with the SPOC PRF and its validation"),
    ("05-transit-fits-and-false-positive-probabilities", "Transit fits and false-positive probabilities",
     "MCMC, NASA's pipeline history, TRICERATOPS"),
    ("06-candidates", "Candidates", "four candidates and one false positive in detail"),
    ("07-weak-signals", "Weak signals", "21 weak signals and a pixel re-check"),
    ("08-limitations-and-reproducibility", "Limitations and reproducibility", "caveats, commands, compute, data"),
]


def page(title, body, rel="", active="", math=False, description=DESCRIPTION):
    nav = [("Overview", f"{rel}index.html", "overview"), ("Write-up", f"{rel}writeup/index.html", "writeup"),
           ("Follow-up", f"{rel}follow-up.html", "follow-up"), ("Code", REPO, "code")]
    links = "".join(f'<a href="{h}" class="{"on" if k == active else ""}">{t}</a>' for t, h, k in nav)
    mathjax = ""
    if math:
        mathjax = ('<script>window.MathJax={tex:{inlineMath:[["\\\\(","\\\\)"]],displayMath:[["\\\\[","\\\\]"]]},'
                   'options:{processHtmlClass:"arithmatex"}};</script>'
                   '<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>')
    full_title = title if title == TITLE else f"{title} · TESS M-dwarf deep search"
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(full_title)}</title>
<meta name="description" content="{html.escape(description)}">
<meta property="og:title" content="{html.escape(full_title)}">
<meta property="og:description" content="{html.escape(description)}">
<meta property="og:image" content="{SITE}/assets/figures/fig07_transits.png">
<meta property="og:type" content="website">
<link rel="icon" href="{rel}assets/favicon.svg" type="image/svg+xml">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&family=JetBrains+Mono:wght@400&display=swap" rel="stylesheet">
<link rel="stylesheet" href="{rel}assets/style.css">
{mathjax}
</head>
<body>
<header class="bar"><div class="inner"><a class="brand" href="{rel}index.html">TESS M-dwarf deep search</a><nav>{links}</nav></div></header>
{body}
<footer><div class="wide"><span>TESS M-dwarf deep search · {AUTHOR} · 2026</span>
<span><a href="{REPO}">Source on GitHub</a> · MIT License</span></div></footer>
</body>
</html>
"""


# ------------------------------------------------------------------ markdown chapters


def gh(path):
    """GitHub URL for a repository path ('/' suffix = directory)."""
    kind = "tree" if path.endswith("/") or "." not in Path(path).name else "blob"
    return f"{REPO}/{kind}/{BRANCH}/{path}"


def rewrite_links(htm, src_dir, rel):
    """Make Markdown links work on the site: chapters -> .html, figures -> assets, anything else -> GitHub."""
    def fix(m):
        attr, url = m.group(1), m.group(2)
        if url.startswith(("http", "#", "mailto:")):
            return m.group(0)
        anchor = ""
        if "#" in url:
            url, anchor = url.split("#", 1)
            anchor = "#" + anchor
        target = (src_dir / url).resolve()
        try:
            relp = target.relative_to(ROOT).as_posix() + ("/" if url.endswith("/") else "")
        except ValueError:
            return m.group(0)
        if relp.startswith("docs/figures/"):
            new = f"{rel}assets/figures/{Path(relp).name}"
        elif relp.startswith("docs/") and relp.endswith(".md"):
            stem = Path(relp).stem
            new = f"{rel}writeup/{'index' if stem == 'README' else stem}.html"
        elif relp == "FOLLOW_UP.md":
            new = f"{rel}follow-up.html"
        else:
            new = gh(relp)
        return f'{attr}="{new}{anchor}"'
    return re.sub(r'(href|src)="([^"]+)"', fix, htm)


def md_to_html(text):
    return markdown.markdown(text, extensions=["tables", "fenced_code", "toc", "attr_list", "pymdownx.arithmatex"],
                             extension_configs={"pymdownx.arithmatex": {"generic": True},
                                                "toc": {"permalink": False}})


def build_chapters():
    (OUT / "writeup").mkdir(parents=True, exist_ok=True)
    for i, (stem, title, _) in enumerate(CHAPTERS):
        text = (ROOT / "docs" / f"{stem}.md").read_text()
        body = rewrite_links(md_to_html(text), ROOT / "docs", "../")
        body = re.sub(r"<table>", '<div class="tablewrap"><table>', body).replace("</table>", "</table></div>")
        prev_ = CHAPTERS[i - 1] if i > 0 else None
        next_ = CHAPTERS[i + 1] if i + 1 < len(CHAPTERS) else None
        pager = '<nav class="pager">'
        pager += (f'<a href="{prev_[0]}.html"><span class="k">Previous</span>{prev_[1]}</a>' if prev_ else "<span></span>")
        pager += (f'<a href="{next_[0]}.html" style="text-align:right"><span class="k">Next</span>{next_[1]}</a>'
                  if next_ else '<a href="index.html" style="text-align:right"><span class="k">Back to</span>Contents</a>')
        pager += "</nav>"
        (OUT / "writeup" / f"{stem}.html").write_text(
            page(title, f'<main class="wrap"><article>{body}{pager}</article></main>', rel="../", active="writeup",
                 math=True))
    items = "".join(f'<li><a href="{s}.html">{i + 1}. {t}</a><span>{d}</span></li>'
                    for i, (s, t, d) in enumerate(CHAPTERS))
    body = (f'<main class="wrap"><article><h1>Technical write-up</h1>'
            f'<p class="lead">A complete account of the search, its calibration and its results, in eight chapters. '
            f'Each stands on its own; together they document every step from raw TESS data to the four candidates.</p>'
            f'<ol class="chapters">{items}</ol>'
            f'<p class="small muted" style="margin-top:22px">Also available as Markdown in the repository: '
            f'<a href="{gh("docs/")}">docs/</a>. Per-candidate dossiers: <a href="{gh("results/candidates/")}">'
            f'results/candidates/</a>. Generated summary of every number: <a href="{gh("REPORT.md")}">REPORT.md</a>.</p>'
            f'</article></main>')
    (OUT / "writeup" / "index.html").write_text(page("Technical write-up", body, rel="../", active="writeup"))


# ------------------------------------------------------------------ overview page


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def fig(name, caption, alt):
    return (f'<figure><div class="card"><img src="assets/figures/{name}" alt="{html.escape(alt)}" loading="lazy">'
            f'</div><figcaption>{caption}</figcaption></figure>')


def build_index():
    from tess_search.assessments import facts

    summ = json.loads((ROOT / "results" / "summary.json").read_text())
    cs = read_csv(ROOT / "results" / "candidates" / "summary.csv")
    loc = read_csv(ROOT / "results" / "hardening" / "localization_summary.csv")
    planets = [r for r in loc if r["kind"] == "planet"]
    nebs = [r for r in loc if r["kind"] == "neb"]
    inj = [r for r in loc if r["kind"] == "injection" and r["reliable"] == "True"]
    inj_ok = sum(r["recovered_correct"] == "True" for r in inj)

    rows = []
    for r in cs:
        tic = int(r["tic"])
        cand = r["verdict"].startswith("Candidate")
        name = 'TOI-218 <span class="muted">· TIC 32090583</span>' if tic == 32090583 else f"TIC {tic}"
        period = float(r["period"])
        per = f"{period:.4f} d" + (f" ({period * 24:.1f} h)" if period < 1 else "")
        rp = f'{float(r["rp_rearth"]):.2f} ± {float(r["rp_err"]):.2f}' if cand else "–"
        teq = f'{float(r["teq_k"]):.0f} K' if cand else "–"
        fpp = f'{float(r["fpp_cleared"]):.2f}' if cand and r["fpp_cleared"] else "–"
        status = {32090583: "new signal", 229689348: "NASA TCE", 149390648: "NASA TCE", 198412174: "NASA TCE",
                  294053492: "eclipsing binary"}[tic]
        rows.append(f'<tr class="{"" if cand else "fp"}"><td>{name}</td><td class="num">{per}</td>'
                    f'<td class="num">{rp}</td><td class="num">{teq}</td>'
                    f'<td class="num">{float(r["loc_offset_arcsec"]):.1f} ± {float(r["loc_err_arcsec"]):.1f}″</td>'
                    f'<td class="num">{fpp}</td><td>{status}</td></tr>')
    table = ('<div class="tablewrap"><table><thead><tr><th>Target</th><th>Period</th><th>Radius (R⊕)</th>'
             '<th>T<sub>eq</sub></th><th>Source offset</th><th>FPP</th><th>Status</th></tr></thead><tbody>'
             + "".join(rows) + "</tbody></table></div>")

    f218, f229, f149, f198, f294 = (facts(t) for t in (32090583, 229689348, 149390648, 198412174, 294053492))
    cards = f"""
<div class="cards">
 <div class="cand"><h3>TOI-218, third signal <span class="tag">strongest</span></h3>
  <div class="meta">P = {f218['period']:.4f} d · R<sub>p</sub> = {f218['rp']} R⊕ · T<sub>eq</sub> ≈ {f218['teq']:.0f} K · FPP {f218['fpp_c']:.3f}</div>
  <ul><li>A third periodic signal beside TOI-218.01 and .02, in 437 transits and both halves of the data.</li>
  <li>TOI-218 has an equal-brightness wide-binary companion 13.5″ away; the pixels place all three signals on TOI-218 and exclude the companion at {f218['comp_sig']:.1f}σ for the new one.</li>
  <li>With the existing Gemini speckle imaging, FPP {f218['fpp_cc']}, below TRICERATOPS's validation thresholds.</li>
  <li>Open: flaring host; no ground-based light curve has yet seen the transit.</li></ul></div>
 <div class="cand"><h3>TIC 229689348 <span class="tag">11.2-hour orbit</span></h3>
  <div class="meta">P = {f229['period']:.4f} d · R<sub>p</sub> = {f229['rp']} R⊕ · T<sub>eq</sub> ≈ {f229['teq']:.0f} K · FPP {f229['fpp_c']:.3f}</div>
  <ul><li>Flat-bottomed, {f229['n_transits']:,} transits; flagged three times by NASA's pipeline but never promoted.</li>
  <li>NASA's 55″ source offset came from difference images that failed its own quality metric; the joint localization puts the source on the target ({f229['loc_off']:.1f} ± {f229['loc_err']:.1f}″).</li>
  <li>Open: unresolved bound companions, which imaging would test.</li></ul></div>
 <div class="cand"><h3>TIC 149390648</h3>
  <div class="meta">P = {f149['period']:.4f} d · R<sub>p</sub> = {f149['rp']} R⊕ · T<sub>eq</sub> ≈ {f149['teq']:.0f} K · FPP {f149['fpp_c']:.3f}</div>
  <ul><li>Earth-sized; on target in a crowded field, in agreement with NASA's own offset.</li>
  <li>NFPP {f149['nfpp_c']:.5f} once the 15 neighbours the pixels exclude are cleared.</li>
  <li>Open: the first half of the data alone does not lock onto the period.</li></ul></div>
 <div class="cand"><h3>TIC 198412174 <span class="tag">host ambiguous</span></h3>
  <div class="meta">P = {f198['period']:.4f} d · R<sub>p</sub> = {f198['rp']} R⊕ · T<sub>eq</sub> ≈ {f198['teq']:.0f} K · FPP {f198['fpp_c']:.3f}</div>
  <ul><li>Near-grazing (b ≈ {f198['b']:.2f}); on or near the target.</li>
  <li>A T = 18.2 star {f198['near_sep']:.1f}″ away is below the localization's resolution, and TRICERATOPS weighs a planet on an unseen companion ({f198['stp']:.0f}%) about equally with the target ({f198['tp']:.0f}%).</li>
  <li>Open: needs adaptive-optics or speckle imaging.</li></ul></div>
</div>"""

    stats = (f'<div class="stats"><div class="stat"><div class="v">{summ["stars"]:,}</div><div class="k">M dwarfs, 20–44 sectors each</div></div>'
             f'<div class="stat"><div class="v">{summ["signals_total"]:,}</div><div class="k">periodic signals detected</div></div>'
             f'<div class="stat"><div class="v">25 / 26</div><div class="k">confirmed transiting planets recovered</div></div>'
             f'<div class="stat"><div class="v">4</div><div class="k">Earth-sized candidates in no catalogue</div></div></div>')

    items = "".join(f'<li><a href="writeup/{s}.html">{i + 1}. {t}</a><span>{d}</span></li>' for i, (s, t, d) in enumerate(CHAPTERS))
    bib = f"""@misc{{tess_mdwarf_deep_search_2026,
  author       = {{{AUTHOR}}},
  title        = {{{TITLE}}},
  year         = {{2026}},
  howpublished = {{\\url{{{REPO}}}}},
  note         = {{Version of October 2026}}
}}"""

    body = f"""
<main>
<div class="wide hero">
  <div class="eyebrow">TESS · M dwarfs · transit search</div>
  <h1 class="title">{TITLE}</h1>
  <p class="subtitle">{DESCRIPTION}</p>
  <div class="byline"><span><a href="{AUTHOR_URL}">{AUTHOR}</a></span><span>October 2026</span></div>
  <div class="buttons">
    <a class="btn" href="{REPO}">{ICON['code']}Code</a>
    <a class="btn" href="writeup/index.html">{ICON['doc']}Technical write-up</a>
    <a class="btn" href="{gh('results/candidates/')}">{ICON['list']}Candidate dossiers</a>
    <a class="btn" href="follow-up.html">{ICON['scope']}Follow-up</a>
    <a class="btn" href="{gh('results/followup_planning/transits.csv')}">{ICON['table']}Transit predictions</a>
  </div>
</div>

<div class="wide">{fig('fig07_transits.png', '<b>The four candidates.</b> Phase-folded TESS photometry in 8-minute bins with the median MCMC transit model. Radii include the stellar-radius uncertainty.', 'Folded transit light curves of the four candidates')}</div>

<section class="wrap" id="summary">
  <h2>Summary</h2>
  <p class="lead">M dwarfs near the TESS continuous viewing zones now have years of 2-minute photometry, enough to
  reach Earth-sized planets. This work searched the {summ['stars']:,} M dwarfs with at least 20 sectors of SPOC
  light curves through Sector 107 and found four Earth-sized transit candidates that appear in no planet or
  candidate catalogue: a third signal in the TOI-218 system and three signals NASA's pipeline had flagged but
  never promoted, including an 11.2-hour orbit around TIC 229689348.</p>
  <p>Each candidate was then tested beyond the light curve. A new pixel-level localization, validated on confirmed
  planets, known nearby eclipsing binaries and eclipses planted in the real images, places all four on their
  target stars and shows that a fifth first-pass candidate is an eclipsing binary 22″ away. With the neighbours
  the pixels exclude treated as cleared, all four meet the TRICERATOPS <i>likely planet</i> criteria, and TOI-218's
  new signal, the only one with existing high-resolution imaging, also falls below the validation thresholds. None is
  yet confirmed; the remaining questions are exactly what ground-based photometry and high-resolution imaging
  answer.</p>
  {stats}
</section>

<section class="wide" id="candidates">
  <h2>Candidates</h2>
  {table}
  <p class="small muted">Source offset: where the light goes missing relative to the target, from the joint
  pixel-level fit (1.5″ systematic floor included). FPP: TRICERATOPS false-positive probability with the neighbours
  the pixels exclude treated as cleared, before high-resolution imaging (TOI-218 with its existing imaging:
  FPP {f218['fpp_cc']}). NASA TCE: flagged by NASA's SPOC pipeline
  as a Threshold Crossing Event but never promoted to a TESS Object of Interest.</p>
  {cards}
</section>

<section class="wrap" id="search">
  <h2>Search and vetting</h2>
  <p>Light curves are cleaned of flares, detrended with a robust biweight filter, and searched from 0.4 to 40 days.
  Because eight years of data with year-long gaps would need a prohibitively fine coherent period grid, box least
  squares runs on each observing season separately and the seasons' likelihoods are added: a real planet adds up
  in every season, noise does not. Eighteen tests then look for eclipsing binaries, contamination, artefacts and
  stellar variability, with noise measured at the transit's own timescale.</p>
  {fig('fig02_funnel.png', '<b>From signals to candidates.</b> Signals remaining after each stage.', 'Bar chart of signals remaining after each stage')}
  <p>The pipeline recovers 25 of 26 confirmed transiting planets around these stars, keeps 74% of 300 planets
  injected into the real light curves (89% of 2–4 R⊕ planets inside 15 days), and produces no false candidates from
  200 light curves turned upside down.</p>
  {fig('fig03_completeness.png', '<b>Completeness.</b> Fraction of injected planets found and kept, by radius and period.', 'Completeness as a function of planet radius')}
</section>

<section class="wrap" id="localization">
  <h2>Where the light goes missing</h2>
  <p>TESS pixels are 21″ wide, so an eclipsing binary a few pixels away can leak a planet-sized dip into the
  target's light curve. For every sector the images taken during transit are subtracted from those just before
  and after, and all sectors are fitted together with NASA's pixel response function, calibrated per sector on Gaia
  stars, to find where the light disappeared. Per-pixel errors come from fake transits at random times; position
  errors include a 1.5″ floor measured on {len(inj) + len(planets)} sources of known position.</p>
</section>
<div class="wide">{fig('fig06_localization_maps.png', '<b>Localization of the five signals.</b> Shading shows how strongly each position is excluded as the source (dark: allowed); the dashed contour is the 3σ region; circles are Gaia DR3 stars. TIC 294053492’s light loss lies on a faint star 22″ north-east.', 'Localization maps for five signals')}</div>
<section class="wrap" style="border-top:0;padding-top:0">
  <p>The method was validated before it was trusted: all {len(planets)} confirmed planets tested come out on their
  own star, all {len(nebs)} signals the TESS Follow-up Observing Program had traced to nearby eclipsing binaries come
  out off target, and {inj_ok} of {len(inj)} reliable synthetic eclipses planted in the real pixels are traced to the
  correct star.</p>
  {fig('fig04_localization_validation.png', '<b>Validation.</b> Left: position errors of sources with known positions against the expectation. Right: how strongly the target is excluded, for confirmed planets, known nearby eclipsing binaries and this work’s signals.', 'Localization validation')}
</section>

<section class="wrap" id="pipeline">
  <h2>What NASA's pipeline saw</h2>
  <p>Three candidates were SPOC Threshold Crossing Events that never became TOIs. Their reported SNR fell as data
  accumulated. Folding this work's light curves at each SPOC period reproduces the drop: a period error of a few
  10<sup>−5</sup> d smears a sub-hour transit by hours over a 2,000-day baseline. The signals did not fade.</p>
  {fig('fig08_spoc_period_drift.png', '<b>Period drift.</b> SNR of this work’s light curve versus trial period, with the periods SPOC adopted in each run.', 'SNR versus trial period')}
</section>

<section class="wrap" id="statistics">
  <h2>How likely is each alternative?</h2>
  <p>TRICERATOPS weighs a planet on the target against eclipsing binaries, unresolved companions, background stars
  and resolved neighbours. It has no pixel information, so it is run twice: on TESS photometry alone, and with the
  neighbours the localization excludes treated as cleared. For TOI-218 this removes the twin companion as a host
  and lowers the FPP from {f218['fpp']:.2f} to {f218['fpp_c']:.2f}. What remains for every candidate is unresolved bound
  companions, which high-resolution imaging tests directly. TOI-218 already has such imaging: Gemini-South speckle
  observations from 2020, taken for its two known TOIs. With that contrast curve, FPP {f218['fpp_cc']}
  and NFPP {f218['nfpp_cc']}, below the thresholds TRICERATOPS uses for validation (FPP &lt; 0.015,
  NFPP &lt; 0.001). It is not called validated here, because the host flares and no ground-based light curve has yet
  seen the transit.</p>
  {fig('fig09_triceratops.png', '<b>Scenario probabilities</b> from TESS photometry alone, with pixel-excluded neighbours cleared, and for TOI-218 with its existing speckle imaging.', 'TRICERATOPS scenario probabilities')}
</section>

<section class="wrap" id="false-positive">
  <h2>A false positive caught</h2>
  <p>TIC 294053492 ({f294['period']:.4f} d) passed every light-curve test and was a candidate in the first pass. In
  the pixels its light loss lies {f294['loc_off']:.1f} ± {f294['loc_err']:.1f}″ north-east of the target, which is
  excluded at {f294['target_sig']:.1f}σ, independently in odd and even sectors. The best-fit position is
  {f294['src_dbest']:.1f}″ from a G = {f294['src_g']:.1f} background star that would need a ~{f294['src_depth']:.0f}%
  eclipse: an ordinary eclipsing binary.</p>
</section>

<section class="wrap" id="weak">
  <h2>Weak signals</h2>
  <p>Twenty-one weaker signals passed the rejecting tests with flags. The inverted light curves predict about 38
  such false alarms in this sample, many cluster at 36–41-day periods where few transits exist, and six are not
  reproduced in the pixels at all.</p>
  {fig('fig10_weak_periods.png', '<b>Weak-signal periods</b> in real data and in inverted data, which contain no planets.', 'Weak signal periods')}
</section>

<section class="wrap" id="follow-up">
  <h2>Follow-up</h2>
  <p>Each candidate's remaining question has a specific observation that answers it: seeing-limited photometry that
  resolves the neighbours during a predicted transit, and adaptive-optics or speckle imaging for unresolved
  companions. Ephemerides, upcoming observable transits from the Las Cumbres Observatory sites, and what each
  observation would settle are on the <a href="follow-up.html">follow-up page</a>.</p>
</section>

<section class="wrap" id="writeup">
  <h2>Technical write-up</h2>
  <ol class="chapters">{items}</ol>
</section>

<section class="wrap" id="reproduce">
  <h2>Reproducing</h2>
  <p>Everything is generated by numbered scripts from public data (about 56 GB from MAST), on a fanless MacBook
  Air under a thermal guard. Commands, seeds and compute times are in
  <a href="writeup/08-limitations-and-reproducibility.html">Chapter 8</a>.</p>
<pre><code>uv venv --python 3.12 .venv &amp;&amp; uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python -m pytest tests/</code></pre>
</section>

<section class="wrap" id="cite">
  <h2>Citation</h2>
<pre><code>{html.escape(bib)}</code></pre>
  <p class="small muted">This work uses TESS data from MAST and SPOC, the TESS Input Catalog, Gaia DR3 via VizieR,
  ExoFOP and the NASA Exoplanet Archive, and the open-source packages Astropy, Astroquery, NumPy, SciPy, pandas,
  Matplotlib, wotan, batman, emcee, lightkurve, TESS_PRF and TRICERATOPS.</p>
</section>
</main>"""
    (OUT / "index.html").write_text(page(TITLE, body, active="overview"))


# ------------------------------------------------------------------ follow-up page


def build_followup():
    fu = ROOT / "FOLLOW_UP.md"
    content = rewrite_links(md_to_html(fu.read_text()), ROOT, "") if fu.exists() else "<h1>Follow-up</h1>"
    rows = read_csv(ROOT / "results" / "followup_planning" / "transits.csv")
    order = ["TIC 32090583 (TOI-218)", "TIC 229689348", "TIC 149390648", "TIC 198412174"]
    by = {}
    for r in rows:
        by.setdefault(r["candidate"], []).append(r)
    parts = ['<h3 id="next-transits">Next observable transits</h3><p>The next eight per candidate (UTC). 1σ is the '
             'mid-transit uncertainty; the window includes 30 minutes of baseline on each side.</p>']
    for cand in sorted(by, key=lambda c: order.index(c) if c in order else 99):
        rs = sorted(by[cand], key=lambda r: r["mid_transit_utc"])[:8]
        body = "".join(f'<tr><td class="num">{r["mid_transit_utc"]}</td><td class="num">±{float(r["sigma_min"]):.0f} min</td>'
                       f'<td class="num">{r["window_start_utc"][11:]}–{r["window_end_utc"][11:]}</td>'
                       f'<td>{html.escape(r["site"].replace("LCO ", ""))}</td><td class="num">{float(r["min_altitude_deg"]):.0f}°</td></tr>'
                       for r in rs)
        parts.append(f'<h4>{html.escape(cand)}</h4>'
                     f'<table><thead><tr><th>Mid-transit</th><th>1σ</th>'
                     f'<th>Window</th><th>Site</th><th>Min. alt.</th></tr></thead><tbody>{body}</tbody></table>')
    tables = "".join(parts)
    marker = "<!-- transit-tables -->"
    content = content.replace(marker, tables) if marker in content else content + tables
    content = re.sub(r"<table>", '<div class="tablewrap"><table>', content)
    content = content.replace("</table>", "</table></div>")
    (OUT / "follow-up.html").write_text(page("Follow-up", f'<main class="wrap"><article>{content}</article></main>',
                                             active="follow-up", math=True))


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "assets" / "figures").mkdir(parents=True)
    shutil.copy(ROOT / "site" / "style.css", OUT / "assets" / "style.css")
    (OUT / "assets" / "favicon.svg").write_text(FAVICON)
    for p in (ROOT / "docs" / "figures").glob("*.png"):
        shutil.copy(p, OUT / "assets" / "figures" / p.name)
    (OUT / ".nojekyll").write_text("")
    build_index()
    build_chapters()
    build_followup()
    print(f"built {sum(1 for _ in OUT.rglob('*.html'))} pages -> {OUT}")


if __name__ == "__main__":
    main()
