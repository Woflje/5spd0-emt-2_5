import numpy as np
parameters = {
    "wavelength": 100,
    "input_mesh_path": "examples/sphere_z_offset.dat",
    "E_field_in": {
        # E-Field Amplitude, Direction [phi,theta], Polarization angle relative to phi and theta
        "amplitude": 1,
        "direction": [np.pi/4,np.pi/4],
        "polarization": np.pi/4
    },
    "E_farfield": {
         # Far-Field Amplitude, Direction [phi,theta], Polarization angle relative to phi and theta
        "amplitude": 1,
        "direction": [np.array([np.pi*2.]), np.linspace(0.001,2*np.pi,180)],
        "polarization": np.array([0])
    },
    "order_dunavant": 5,
    "order_duffy": 5
}