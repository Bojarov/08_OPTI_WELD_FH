
def f_b_field(y, I, z):
    """
    y-component of the magnetic field of infinite straight wire.
    """
    return 0.0000002 * I * z / (y * y + z * z)
