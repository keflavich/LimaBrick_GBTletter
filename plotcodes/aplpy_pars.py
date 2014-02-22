# Default aplpy parameters

def set_aplpy_pars(F):
    F.set_auto_refresh(False)
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    F.recenter(0.26, 0.035, height=0.2, width=0.2)
    F.refresh()
