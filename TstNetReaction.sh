#! /bin/sh
python net_balance.py  ../inout/camx403.20000830.base5b.psito2n2.1km.clinton.1km.hourly.cycle.8.18.ext \
  ../inout/Net_balance_clinton.txt

python net_balance.py  ../inout/camx403.20000830.base5b.psito2n2.1km.source.1km.hourly.cycle.8.18.ext \
  ../inout/Net_balance_source.txt
