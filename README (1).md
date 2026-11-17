[![make.fun · gta6-twitch-ccu-forecast-01](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/Hawdiias/gta6-twitch-forecast/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/Hawdiias/gta6-twitch-forecast/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)

# GTA VI Launch-Day Twitch Forecaster

A forecasting program for the make.fun bounty "Forecast GTA VI's launch-day
Twitch viewers hour by hour". Given a case file holding everything recorded
up to the end of hour 1 of a game's release day — that hour's viewer count
plus the hourly curves of dozens of past big launches — it writes, for each
of hours 2 to 24, five quantiles (5%, 25%, 50%, 75%, 95%) of that hour's
concurrent Twitch viewers in the game's category.

## Method

Every past launch in the case contributes its release-day log-ratio curve,
ln(viewers in hour h / viewers in hour 1). A launch's vote is weighted by
how close its own hour-1 audience was to this case's hour-1 audience — a
Gaussian kernel on the log10 scale with a half-dex bandwidth — because
launches of a similar size trace the shape a new launch is likely to
follow: the slow build before the unlock wave, the peak, the long decay.
The forecast for hour h is the case's hour-1 count times the exponential
of the weighted quantile of those log-ratios, with the far tails trimmed
at 2.2 log units from the median so one freak curve cannot stretch the
spread. Launches whose release day followed a paid early-access window are
left out of the fit: their day-one shape is not a launch shape. Hours with
too little weighted evidence fall back to the unweighted pool, then to a
persistence forecast with a widening spread.

The program is a pure function of its case file: no clock, no network, no
randomness — two runs write byte-identical forecasts.

## Running

The harness runs it as `python3 main.py` (Python 3.13, standard library
only). It reads every `*.in` case in `$MAKEFUN_INPUT_DIR` (else `/input`)
and writes one `<id>.out` per case to `$MAKEFUN_OUTPUT_DIR` (else
`/output`).

## License

MIT — see LICENSE.
