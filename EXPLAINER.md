# Plain-language summary

## In one sentence
This project searches eight years of NASA telescope data on 1,279 small, cool stars for the tiny, regular
dimming caused by planets passing in front of them, and then checks every signal it finds, in the data
and in the raw images, to weed out the many things that are not planets.

## The telescope
**TESS** (Transiting Exoplanet Survey Satellite) is a NASA space telescope launched in 2018. It watches a
patch of sky for about a month (a "sector"), measuring the brightness of tens of thousands of stars every
2 minutes, then moves on. Some parts of the sky, near the poles of Earth's orbit, are watched again and
again. By Sector 107 some stars have been watched for over 40 sectors, more than three years of total
observing.

## How to find a planet you cannot see
When a planet passes in front of its star (a **transit**), it blocks a little of the star's light. The
star dims slightly for an hour or two, then recovers, and the dip repeats every orbit. An Earth-sized
planet in front of the Sun dims it by 0.008%; in front of a small **M dwarf** (a red dwarf, a third of the
Sun's size) the same planet dims it by about 0.1%, ten times easier to see. That is why this search
targets M dwarfs.

## The steps
1. **Pick stars.** From 1.7 million TESS light-curve files, the 1,279 bright M dwarfs that TESS watched for
   at least 20 sectors.
2. **Clean the data.** Stars flicker (starspots rotate in and out of view) and flare (red dwarfs throw off
   bright flashes). Flares are removed and the slow wiggles flattened, so only short dips remain.
3. **Search.** For every orbital period from 9.6 hours to 40 days, the search checks whether the brightness
   dips at regular intervals ("box least squares"). Eight years of data with big gaps is a lot of
   combinations, so each observing year is searched separately and the results are added up: a real
   planet shows up in every year, random noise does not.
4. **Vet.** Most dips are not planets. Two stars orbiting each other (an eclipsing binary) make dips too,
   usually of alternating depth or with a second, smaller dip halfway between. A neighbouring star can
   leak its own eclipses into the measurement. Spacecraft glitches repeat on TESS's 13.7-day orbit. Every
   signal goes through eighteen such tests.
5. **Is it already known?** Every survivor is compared with NASA's lists of known planets, official TESS
   candidates, community candidates, known eclipsing binaries, and every signal NASA's own pipeline ever
   flagged.

## Evidence that it works
* **It rediscovers known planets.** It finds 25 of the 26 known transiting planets in its range on its
  own, including all four planets of **TOI-700** (two of them Earth-sized and in or near the zone where
  liquid water could exist) and the three transiting planets of **L 98-59**. Its planet sizes agree with
  published ones to within about one and a half standard deviations (TOI-700 d: 1.18 Earth radii here,
  1.07 published; L 98-59 c: 1.34 here, 1.39 published).
* **Fake planets.** 300 synthetic planets were planted into real data: about three in four came back out,
  nearly nine in ten for planets twice Earth's size, and one in five for planets smaller than Earth.
* **Upside-down data.** 200 light curves were flipped upside down (dips become bumps) and searched again.
  Anything found in flipped data is fake by construction, and not one flipped star produced a
  "candidate", so the candidates are unlikely to be noise.

## The second, deeper check
The first pass found five Earth-sized signals that passed every test on the light curves. A light curve,
though, is just one number per moment: the total light in a small box of pixels. TESS pixels are large
(21 arcseconds), so several stars can share that box. If one of the *other* stars is a pair of stars
eclipsing each other, its dips leak into the box and look like a small planet. So the second pass went
back to the raw images.

* **Where does the light go missing?** For every observing month, the images taken during the dips are
  subtracted from the images just before and after, leaving only the light that disappeared. Fitting all
  months together with NASA's model of how a star spreads over the pixels gives the position on the sky
  where the light went missing. Before being trusted, the method was tested on known cases: confirmed
  planets come out on their own star, signals the TESS team had already traced to neighbouring stars come
  out on those neighbours, and fake eclipses planted into the real images were traced back to the right
  star 34 times out of 35.
* **How likely is each alternative?** TRICERATOPS, a tool used by the TESS team, weighs "planet on this
  star" against every way the dip could be faked (a binary star behind it, an unseen companion, a
  neighbour) and gives a probability.
* **What did NASA's own pipeline see?** Three of the signals had been flagged by NASA's pipeline but never
  promoted. Its reports show why its confidence seemed to drop over time: its estimate of the orbital
  period was slightly off, which blurs a short dip when years of data are stacked. At the right period the
  signals are as strong as ever.
* **Better sizes.** A careful fit that includes the uncertainty in each star's own size, checked on two
  well-studied planets.

## What it found
* **One of the five was not a planet.** For TIC 294053492 the light goes missing about 22 arcseconds from
  the star, on a faint background star that is almost certainly an eclipsing binary. This is exactly the
  kind of mistake the second pass exists to catch.
* **Four remain candidates**, all on their own star:
  * **TOI-218, a third signal** (every 2.15 days, about Earth's size). TOI-218 turned out to be one of a
    pair of twin red dwarfs orbiting each other far apart; the image check shows the new signal, and the
    two already-known ones, come from TOI-218 itself and not from its twin. A sharp image of the star from
    an 8-metre telescope, taken in 2020 for the known signals, rules out most remaining alternatives: the
    statistical test then puts the chance of a false positive at about 1 in 70,000.
  * **TIC 229689348**, a dip every **11.2 hours** (1.25 Earth sizes, about 1,100 K). NASA's pipeline
    thought the source might be a star 55 arcseconds away; the image check, done at the correct period
    with all the data, puts it on the target.
  * **TIC 149390648** (Earth-sized, 2.84 days) in a crowded patch of sky.
  * **TIC 198412174** (1.35 Earth sizes, 1.39 days): a faint star only 5 arcseconds away is too close to
    rule out, and the dip's shape also fits a planet around a small unseen companion star. A sharp image
    from a large telescope would settle it.
* All four meet TRICERATOPS's "likely planet" bar. "Statistically validated" needs a very low
  false-positive probability plus sharp images from a large telescope; TOI-218's new signal, the only one
  with such images already, passes that numerical bar, but a telescope on the ground still needs to see
  one of its transits, because its star flares. The other three are documented candidates, ready for
  that step.

## What a "candidate" means, and does not mean
A candidate passed every automated test and is not in any catalogue. It is *not* a confirmed planet.
Confirmation needs more observations: telescopes on the ground to rule out nearby eclipsing stars, sharper
images, and measurements of the star's wobble. What a candidate *is*: a specific, reproducible, documented
signal worth those follow-up observations. NASA's ExoFOP follow-up database lists such signals as
"community TOIs", but since August 2026 only after they are published in a refereed journal
([FOLLOW_UP.md](FOLLOW_UP.md) explains the route).

## Glossary
* **Light curve:** a star's brightness over time.
* **Transit depth:** how much the star dims, in ppm (parts per million; 1,000 ppm = 0.1%).
* **SNR (signal-to-noise ratio):** how much the signal stands out from the random flicker. Planets are
  usually only believed above about 7.
* **Red noise:** flicker that is correlated in time (slow wiggles), which fools simple statistics; the
  vetting measures it directly instead of assuming it away.
* **False positive:** a signal that looks like a planet but is something else.
* **Eclipsing binary:** two stars orbiting each other that take turns blocking each other's light.
* **TOI:** TESS Object of Interest, NASA's official list of planet candidates.
* **Difference image:** an image of only the light that disappeared during a dip.
* **FPP (false-positive probability):** the chance a signal is not a planet on the target star.
