try:
    from guis.Simplewx import StartWx
    StartGUI = StartWx
except:
    from guis.SimpleTk import StartTk
    StartGUI = StartTk
    