from __future__ import print_function
__all__ = ['IRRTable', 'NetTableMaker', 'PhyTableMaker', 'PtbMaker', 'SumTableMaker']

from . import IRRTable
from . import NetTableMaker
from . import PhyTableMaker
from . import PtbMaker
from . import SumTableMaker

def net_balance(mechanism, mrg_data_path, output_dir):
    from permm.analyses.net_balance.PhyTableMaker import PhyTable
    from permm.analyses.net_balance.SumTableMaker import SumTable
    from permm.analyses.net_balance.NetTableMaker import NetTables
    from permm.analyses.net_balance.PhyTableMaker import PhyTables
    from permm.analyses.net_balance.PhyTableMaker import VOCTable
    from permm.analyses.net_balance.PtbMaker import PtbTable
    from permm.analyses.net_balance.IRRTable import IRRTable
    from permm import getmech
    from permm.netcdf import NetCDFFile
    
    get_prepared_mech = getmech.get_prepared_mech
    get_pure_mech = getmech.get_pure_mech
    try:
        mech = get_prepared_mech(mechanism)
    except:
        mech = get_pure_mech(mechanism)
    

    if output_dir is None:
        net_data_path = '.'.join(mrg_data_path.split('.')[:-1])+'.net.'
    else:
        net_data_path = output_dir+'.'
    
    mrg_file = NetCDFFile(mrg_data_path, 'r')
    mech.set_mrg(mrg_file)

    print(SumTable(mech), file=file(net_data_path+'sum','wb'))
    print(NetTables(mech), file=file(net_data_path+'net','wb'))
    print(PhyTables(mech), file=file(net_data_path+'phy','wb'))
    print(VOCTable(mech), file=file(net_data_path+'voc','wb'))
    print(PtbTable(mech), file=file(net_data_path+'ptb','wb'))
    mech = get_pure_mech(mechanism)
    mech.set_mrg(mrg_file)

    
    print(IRRTable(mech), file=file(net_data_path+'irr','wb'))