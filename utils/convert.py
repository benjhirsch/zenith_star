def is_float(x):
    try:
        y = float(x) 
        return True
    except:
        return False

def to_vmag(star, mag_cols, mag_conversion):
    try:
        import simpleeval
    except ImportError:
        print('Cannot convert catalog magnitude values to visual magnitude without simpleeval. Run "pip install simpleeval" from the command line.')
        return
    
    se = simpleeval.SimpleEval()

    for var in mag_cols:
        se.names[var] = star[var]

    if all([is_float(star[var]) for var in mag_cols]):
        vmag = se.eval(mag_conversion)
    elif is_float(star[mag_cols[0]]):
        vmag = star[mag_cols[0]]

    return vmag

def to_bv(star, mag_cols, bv_conversion):
    try:
        import simpleeval
    except ImportError:
        print('Cannot convert catalog magnitude values to B-V color index without simpleeval. Run "pip install simpleeval" from the command line.')
        return
    
    se = simpleeval.SimpleEval()
    
    for var in mag_cols:
        se.names[var] = star[var]

    if all([is_float(star[var]) for var in mag_cols]):
        bv = se.eval(bv_conversion)
    else:
        bv = None

    return bv

def to_angle(theta, theta_format):
    from astropy.coordinates import Angle
    #handles hms/dms strings
    return Angle(theta if isinstance(theta, float) else tuple(float(a) for a in theta.split()), unit=theta_format)

def make_float(x):
    try:
        return float(x)
    except:
        return x
    
def star_id(obj_id, star):
    #format star id string
    star_id_str = obj_id['format']
    for var in obj_id['columns']:
        star_id_str = star_id_str.replace('{%s}' % var, str(star[var]))
    return star_id_str