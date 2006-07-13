#! /bin/sh

## $LastChangedDate: 2006-04-24 11:56:17 -0400 (Mon, 24 Apr 2006) $
## $LastChangedRevision: 115 $
## $LastChangedBy: jeffries $


python net_balance_CB4.py  ../inout/X420U.000830.b1b.TQ.std.clinton.hourly.cycle.8.18.ext \
  ../inout/Clinton_8h_000830_0818h.txt

python net_balance_CB4.py  ../inout/X420U.000825.b1b.TQ.std.bayland.hourly.cycle.8.18.ext \
  ../inout/Bayland_8h_000825_0818h.txt

python net_balance_CB4.py  ../inout/X420U.000825.b1b.TQ.std.HALC.hourly.cycle.8.18.ext \
  ../inout/HALC_8h_000825_0818h.txt

python net_balance_CB4.py  ../inout/X420U.000830.b1b.TQ.std.HALC.hourly.cycle.8.18.ext \
  ../inout/HALC_8h_000830_0818h.txt
  
python net_balance_CB4.py  ../inout/X420U.000830.b1b.TQ.std.DRPK.hourly.cycle.8.18.ext \
  ../inout/DRPK_8h_000830_0818h.txt  
