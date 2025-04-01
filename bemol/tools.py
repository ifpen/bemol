import numpy as np


def calculateVelocity(wind:float,omega:float,rad:float,azi:float,yaw:float,tilt:float,precone:float):
    """Calculate relative velocity for a given wind configuration

    Parameters
    ----------
    wind : float
        wind speed, m/s
    omega : float
        rotation velocity, rad/s
    rad : float
        radius
    azi : float
        azimuthal angle, radians
    yaw : float
        yaw angle, radians
    tilt : float
        tilt angle, radians
    precone: float
        precone angle, radians

    """
    Ux = wind*(
            (np.cos(yaw)*np.sin(tilt)*np.cos(azi)+np.sin(yaw)*np.sin(azi))*np.sin(precone)
            + np.cos(yaw)*np.cos(tilt)*np.cos(precone)
        )
    Uy = wind*(
            np.cos(tilt)*np.sin(precone)*np.sin(azi)-np.sin(yaw)*np.cos(azi)
        ) + omega*rad*np.cos(precone)
    return Ux, Uy