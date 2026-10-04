# What this project does, in plain language

## The one-sentence version
I wrote software that searches 8 years of NASA telescope data on 1,279 small, cool
stars for the tiny, regular dimming caused by planets passing in front of them, and
then automatically checks every signal it finds to weed out the many things that are
not planets.

## The telescope
**TESS** (Transiting Exoplanet Survey Satellite) is a NASA space telescope launched in
2018. It watches a patch of sky for about a month (a "sector"), measuring the brightness
of tens of thousands of stars every 2 minutes, then moves on. Some parts of the sky,
near the poles of Earth's orbit, get watched again and again. By now (sector 107) some
stars have been watched for over 40 sectors, more than 3 years of total observing.

## How you find a planet you cannot see
When a planet passes in front of its star (a **transit**), it blocks a little of the
star's light. The star dims slightly for an hour or two, then recovers, and the dip
repeats every orbit. An Earth-sized planet in front of the Sun dims it by 0.008%; in
front of a small **M dwarf** star (a red dwarf, a third the Sun's size), the same planet
dims it by about 0.1%, ten times easier to see. That is why this search targets M dwarfs.

## The steps
1. **Pick stars.** From 1.7 million TESS light-curve files I found the 1,279 bright M
   dwarfs that TESS watched for at least 20 sectors.
2. **Clean the data.** Stars flicker (starspots rotate in and out of view) and flare
   (red dwarfs throw off bright flashes). The software removes flares and flattens the
   slow wiggles so only short dips remain.
3. **Search.** For every possible orbital period from 9.6 hours to 40 days, it checks
   whether the brightness dips at regular intervals ("box least squares"). Eight years
   of data with big gaps is a lot of combinations, so it searches each observing year
   separately and adds the results up: a real planet shows up in every year, random
   noise does not.
4. **Vet.** Most dips are not planets. Two stars orbiting each other (an eclipsing
   binary) make dips too, but usually of alternating depth, or with a second smaller dip
   halfway between. A neighbouring star can leak its own eclipses into the measurement;
   then the star's image shifts slightly during the dip. Spacecraft glitches repeat on
   TESS's 13.7-day orbit. The software runs about fifteen such tests on every signal.
5. **Is it already known?** Every survivor is compared against NASA's lists of known
   planets, official TESS candidates, community candidates, known eclipsing binaries,
   and every signal NASA's own pipeline ever flagged.

## How I know it works
* **It rediscovers known planets.** On stars with known planets it found every
  transiting planet in its search range on its own, including all four planets of
  **TOI-700** (two of them Earth-sized and in or near the zone where liquid water could
  exist) and the three transiting planets of **L 98-59**. Its vetting passed them while
  rejecting signals the TESS team had already labelled false positives.
* **Fake planets.** I planted hundreds of synthetic planets into real data and counted
  how many came back out; this measures what sizes and periods the search can detect.
* **Upside-down data.** I flipped every light curve upside down (dips become bumps) and
  searched again. Any "planet" found in flipped data is fake by construction, which
  measures how often the pipeline is fooled by noise.

## What a "candidate" means (and does not mean)
A candidate passed every automated test and is not in any catalogue. It is *not* a
confirmed planet. Confirmation needs more observations: telescopes on the ground to rule
out nearby eclipsing stars, sharper images, and measurements of the star's wobble. What
a candidate *is*: a specific, reproducible, documented signal that is worth those
follow-up observations. Anyone can submit such a signal to NASA's ExoFOP follow-up
database as a "community TOI".

## Words you might be asked about
* **Light curve:** a star's brightness over time.
* **Transit depth:** how much the star dims (in ppm, parts per million; 1000 ppm = 0.1%).
* **SNR (signal-to-noise ratio):** how much the signal stands out from the random
  flicker. Planets are usually only believed above about 7.
* **Red noise:** flicker that is correlated in time (slow wiggles), which fools simple
  statistics; the vetting measures it directly instead of assuming it away.
* **False positive:** a signal that looks like a planet but is something else.
* **TOI:** TESS Object of Interest, NASA's official list of planet candidates.
