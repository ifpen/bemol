import numpy as np


def calculateVelocity(wind:float,omega:float,rad:float,azi:float,yaw:float,tilt:float,precone:float):
    """Calculate relative velocity for a given wind configuration

    Parameters
    ----------
    wind : float
        wind speed in m/s
    omega : float
        rotation velocity
    rad : float
        radius
    azi : float
        azimuthal angle in radians
    yaw : float
        yaw angle in radians
    tilt : float
        tilt angle in radians
    precone: float
        precone angle in radians

    """
    Ux = wind*(
            (np.cos(yaw)*np.sin(tilt)*np.cos(azi)+np.sin(yaw)*np.sin(azi))*np.sin(precone)
            + np.cos(yaw)*np.cos(tilt)*np.cos(precone)
        )
    Uy = wind*(
            np.cos(tilt)*np.sin(precone)*np.sin(azi)-np.sin(yaw)*np.cos(azi)
        ) + omega*rad*np.cos(precone)
    return Ux, Uy