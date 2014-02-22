# Default aplpy parameters

def set_aplpy_pars(F, xc=0.253, yc=0.016, height=0.2, width=0.2, xformat='dd.d', yformat='dd.d'):
    F.set_auto_refresh(False)
    F.set_tick_labels_xformat(xformat)
    F.set_tick_labels_yformat(yformat)
    F.recenter(xc, yc, height=height, width=width)
    F.refresh()
