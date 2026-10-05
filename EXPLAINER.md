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
* **It rediscovers known planets.** It found 25 of the 26 known transiting planets in its
  search range on its own, including all four planets of **TOI-700** (two of them
  Earth-sized and in or near the zone where liquid water could exist) and the three
  transiting planets of **L 98-59**. Its size estimates match the published ones (TOI-700 d:
  1.02 Earth radii here, 1.07 published). Its vetting passed real planets while rejecting
  signals the TESS team had already labelled false positives.
* **Fake planets.** I planted 300 synthetic planets into real data and counted how many
  came back out: about three in four overall, nearly nine in ten for planets twice
  Earth's size, and one in five for planets smaller than Earth.
* **Upside-down data.** I flipped 200 light curves upside down (dips become bumps) and
  searched again. Anything found in flipped data is fake by construction. Not one flipped
  star produced a "candidate", so the candidates below are unlikely to be noise.

## The second, deeper check
The first pass found five Earth-sized signals that passed every test on the light curves.
A light curve, though, is just one number per moment: the total light in a small box of
pixels. TESS pixels are big (21 arcseconds), so several stars can share that box. If one of
the *other* stars is a pair of stars eclipsing each other, its dips leak into the box and
look like a small planet. So the second pass went back to the raw images:

* **Where does the light go missing?** For every observing month, the images taken during
  the dips are subtracted from the images just before and after. What is left shows only
  the light that disappeared. Fitting all months together with NASA's model of how a star
  spreads over the pixels gives the position on the sky where the light went missing.
  Before trusting it on the candidates, it was tested on known cases: confirmed planets
  come out on their own star, signals the TESS team had already traced to neighbouring
  stars come out on those neighbours, and fake eclipses planted into the real images on
  nearby stars were traced back to the right star 34 times out of 35.
* **How likely is each alternative?** A tool called TRICERATOPS, used by the TESS team,
  weighs "planet on this star" against every way the dip could be faked (a binary star
  behind it, an unseen companion, a neighbour) and gives a probability.
* **What did NASA's own pipeline see?** Three of the signals had been flagged by NASA's
  pipeline but never promoted. Its reports showed why its confidence seemed to drop over
  time: its estimate of the orbital period was slightly off, which blurs a short dip when
  years of data are stacked. At the right period the signals are as strong as ever.
* **Better sizes.** A careful fit (MCMC) that includes the uncertainty in each star's
  own size, checked on two well-studied planets.

## What it found
* **One of the five was not a planet.** For TIC 294053492 the light goes missing about
  22 arcseconds away from the star, on a faint background star that is almost certainly
  an eclipsing binary. This is exactly the kind of mistake the second pass exists to catch.
* **Four remain candidates**, all on their own star:
  * **TOI-218, a third signal** (every 2.15 days, about Earth's size). TOI-218 turned out
    to be one of a pair of twin red dwarfs orbiting each other far apart; the image check
    shows the new signal, and the two already-known ones, come from TOI-218 itself and not
    from its twin.
  * **TIC 229689348**, a dip every **11.2 hours** (1.25 Earth sizes, about 1,100 K). NASA's
    pipeline thought the source might be a star 55 arcseconds away; the image check, done at
    the correct period with all the data, puts it on the target.
  * **TIC 149390648** (Earth-sized, 2.84 days) in a crowded patch of sky.
  * **TIC 198412174** (1.35 Earth sizes, 1.39 days): a faint star only 4.7 arcseconds away
    is too close to rule out, and the dip's shape also fits a planet around a small unseen
    companion star. A sharp ground-based image would settle it.
* None is "statistically validated" (that needs a very low false-positive probability plus
  sharp images from a large telescope). They are documented candidates, ready for that step.

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
* **Difference image:** an image of only the light that disappeared during a dip.
* **FPP (false-positive probability):** the chance a signal is not a planet on the target star.
* **Eclipsing binary:** two stars orbiting each other that take turns blocking each other's light.
