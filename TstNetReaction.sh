#! /bin/sh
python net_balance.py  ../inout/camx403.20000830.base5b.psito2n2.1km.clinton.1km.hourly.cycle.8.18.ext \
  ../inout/Clinton_1hour_0830_818_403.txt

python net_balance.py  ../inout/camx420_pa.20000830.base1b.psito2n2.4km.clinton.hourly.cycle.8.18.ext \
  ../inout/Clinton_8hour_0830_818_420.txt

python net_balance.py  ../inout/camx420_pa.20000830.base1b.4time_co.4km.clinton.hourly.cycle.8.18.ext \
  ../inout/Clinton_4timeCO_0830_818_420.txt

python net_balance.py  ../inout/camx420_pa.20000830.base1b.quart_co.4km.clinton.hourly.cycle.8.18.ext \
  ../inout/Clinton_quartCO_0830_818_420.txt

python net_balance.py  ../inout/camx420_pa.20000825.base1b.psito2n2.4km.bayland.hourly.cycle.8.18.ext \
  ../inout/Bayland_8hour_0825_818_420.txt

python net_balance.py  ../inout/camx420_pa.20000825.base1b.4time_co.4km.bayland.hourly.cycle.8.18.ext \
  ../inout/Bayland_4timeCO_0825_818_420.txt

python net_balance.py  ../inout/camx420_pa.20000825.base1b.quart_co.4km.bayland.hourly.cycle.8.18.ext \
  ../inout/Bayland_quartCO_0825_818_420.txt


