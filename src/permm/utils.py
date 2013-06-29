__all__ = ['AttrDict']

import unittest

class AttrDict(dict):
    """For easy access to values"""
    def __new__(cls,*args,**kwds):
        result = dict.__new__(cls,*args,**kwds)
        return result
        
    def __init__(self,*args,**kwds):
        dict.__init__(self, *args, **kwds)
        for k,v in self.iteritems():
            if type(v)==dict:
                self[k]=AttrDict(v)
            else:
                self[k]=v
                
    def __setattr__(self, attr, value):
        self[attr] = value
        
    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            raise AttributeError(attr)


class TestAttrDict(unittest.TestCase):
    def setUp(self):
        self.base_dict = {'one': 1, "two": {"two": 2}}
        self.attr_dict = AttrDict(self.base_dict)
        
    def testInitFromDict(self):
        bd = self.base_dict
        ad = self.attr_dict
        for k in self.base_dict:
            self._assert(bd[k] == ad[k])
            