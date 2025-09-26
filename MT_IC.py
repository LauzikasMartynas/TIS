import numpy as np
import h5py
import shutil
from scipy import interpolate
import MyFunctions as ML

# Input file
source = 'F://isothermal_ellipse_6M_MT.hdf5'

# Temporary copy
path = 'F://E//TIS_temp.hdf5'
shutil.copy(source, path)

# Density distribution
def density_TIS(pos_gas, A, B, c_s):
    G_cgs = 6.67e-8
    if np.shape(np.shape(pos_gas))[0]==1:
        r = np.sqrt(np.sum(pos_gas**2)) * unit_length_cgs
        theta = np.arctan2(np.sqrt(pos_gas[0]**2+pos_gas[1]**2), pos_gas[2])
    else:
        r = np.sqrt(np.sum(pos_gas**2, axis=1)) * unit_length_cgs
        theta = np.arctan2(np.sqrt(pos_gas[:, 0]**2+pos_gas[:, 1]**2), pos_gas[:, 2])
    rho = c_s**2 * A**2 * B * (np.tan(theta/2))**A / (2 * np.pi * G_cgs * r**2 * (
        np.sin(theta))**2 * (1 + B * (np.tan(theta/2))**A)**2)
    return rho

# Units
unit_mass_cgs = 1.989e43            # 1e10 M_Sun
unit_length_cgs = 3.086e21          # kpc
unit_velocity_cgs = 20733964.4      # G=1 in internal units
unit_time_cgs = unit_length_cgs / unit_velocity_cgs
unit_density_cgs = unit_mass_cgs / unit_length_cgs**3
unit_spec_energy_cgs = unit_velocity_cgs**2
unit_pressure_cgs = unit_mass_cgs * unit_velocity_cgs**2 / unit_length_cgs**3

# Constants
G_cgs = 6.67e-8
m_p = 1.67e-24
k_B = 1.38e-16
mu = 0.63
gamma = 5/3
s_in_Myr = 3.16e13

# Bugle + disc
smbh_mass = 0.01 * unit_mass_cgs    # 1e8 M_Sun
sigma = (smbh_mass/(3*unit_mass_cgs/100))**(1/4.4) * 200e5  # from M-sigma
#smbh_mass = 0
r_bulge = 2.0 * unit_length_cgs     # bulge cut-off radius
r_min = 0.01 * unit_length_cgs       # Removes particles closer than this
bulge_mass = 2 * sigma**2 * r_bulge / G_cgs # Bulge mass within r_bulge from isothermal sphere
f_g = 0.10                           # Bulge gas fraction
b_part = 1.0
d_part = 0.0
bulge_gas_mass = bulge_mass * f_g   # Gas mass bulge
ism_temp = 1e4
target_particle_count = 4e5         # target number of particles
particle_mass = (bulge_gas_mass)/target_particle_count/unit_mass_cgs
turb_on = 1                         # turn turbulence on or off
k_large = 6                         # wavenumber limit for turbulence
k_small = 40                        # wavenumber limit for turbulence

f_b = (1.0+f_g)/1.0                        # manual adjust
f_d = (1.0+f_g)/1.0                         # manual adjust
f_v = 1.0                         # manual adjust

seedas = 3
ratio_turb = 0.5

# Gas velocity dispersionsame as stellar one
v_0 =  0.42 * sigma # (np.sqrt(2) - 1) * sigma
c_s = 0.58 * sigma # (2 - np.sqrt(2)) * sigma
#v_0 = 0
#c_s = sigma

# Disc component. v_0_d = 0 disables it.
base = 20
v_0_d = 0#(((base)/(base-1) - np.sqrt(base)/(base-1))) * sigma
c_s_d = 1#(np.sqrt(base)/(base-1) - 1/(base-1)) * sigma


# Gas temperature
gas_spec_energy = c_s**2 / (gamma-1) / unit_spec_energy_cgs * f_b
A = 2 + v_0**2/c_s**2
B = 1

gas_spec_energy_d = c_s_d**2 / (gamma-1) / unit_spec_energy_cgs * f_d
A_d = 2 + v_0_d**2/c_s_d**2
B_d = 1

print("Bulge: A="+str(A)+";B="+str(B))
print("Disc:  A="+str(A_d)+";B="+str(B_d))
print("v_0="+str(v_0/1e5)+"; c_s="+str(c_s/1e5))
print("v_0_d="+str(v_0_d/1e5)+"; c_s_d="+str(c_s_d/1e5))

with h5py.File(path, 'r+') as file:
    header = file['Header']

    # Duomenys
    gas = file['PartType0']

    pos_gas = (gas['Coordinates'][()] - np.array((0.5, 0.5, 0.5)))
    id_gas = gas['ParticleIDs'][()]


    # Sort and filter by radius
    r = np.sqrt(np.sum(pos_gas**2, axis=1)) * unit_length_cgs
    idx_sort = r.argsort()
    pos_gas = pos_gas[idx_sort]
    pos_gas = pos_gas[0:int(target_particle_count), :]
    r = r[idx_sort]
    r = r[0:int(target_particle_count)]
    id_gas = id_gas[idx_sort]
    id_gas = id_gas[0:int(target_particle_count)]
    pos_gas = pos_gas / np.max(r) * r_bulge
    r = np.sqrt(np.sum(pos_gas**2, axis=1)) * unit_length_cgs

    pos_gas = pos_gas[r > r_min]
    r = r[r > r_min]

    rho = b_part * density_TIS(pos_gas, A, B, c_s) * f_g * f_b
    rho_d = d_part * density_TIS(pos_gas, A_d, B_d, c_s_d) * f_g * f_d

    if len(pos_gas) < 0.99 * target_particle_count:
        raise Exception("Source too small.")

    header.attrs['MassTable'] = np.array((0, 0, 0, 0, 0, 0))

    pos = np.vstack((pos_gas))

    # Add rotation
    vel_gas = np.zeros_like(pos_gas)
    vel_gas[:, 0] = pos_gas[:, 1]
    vel_gas[:, 1] = -pos_gas[:, 0]
    vel_norm = np.sqrt(np.sum(vel_gas**2, axis=1))[:, np.newaxis]
    vel_gas = vel_gas / vel_norm * v_0 / unit_velocity_cgs

    # Add rotation for disc
    vel_gas_d = np.zeros_like(pos_gas)
    vel_gas_d[:, 0] = pos_gas[:, 1]
    vel_gas_d[:, 1] = -pos_gas[:, 0]
    vel_norm_d = np.sqrt(np.sum(vel_gas_d**2, axis=1))[:, np.newaxis]
    vel_gas_d = vel_gas_d / vel_norm_d * v_0_d / unit_velocity_cgs * f_v

    vel_gas = vel_gas * (rho/(rho+rho_d))[:, np.newaxis] + vel_gas_d * (rho_d/(rho+rho_d))[:, np.newaxis]

    mass_gas = np.ones_like(pos[:,0]) * particle_mass

    # Determined by the shape of the potential
    T_cs = c_s**2 * mu * m_p / (k_B) * f_b
    T_cs_d = c_s_d**2 * mu * m_p / (k_B) * f_d

    sigma_turb = np.sqrt(T_cs * k_B / (mu * m_p)) / unit_velocity_cgs
    sigma_turb_d = np.sqrt(T_cs_d * k_B / (mu * m_p)) / unit_velocity_cgs

    # Turbulate
    if (np.mean(sigma_turb*turb_on) > 0):
        vel_turb = ML.TurbFieldML(res=256, kmin=k_large*np.pi, method='proj', ratio=ratio_turb, seed=seedas)

        print(f'Interpolating {len(pos)} points.')
        dim = np.linspace(-r_bulge/unit_length_cgs, r_bulge/unit_length_cgs, vel_turb.shape[-1])
        grid = (dim, dim, dim)  # Make position grid same size as turbulent vel field

        # Do linear turbulent velocity interpolation from grid to pos
        vel_g = np.zeros_like(pos_gas)
        for _ in range(3):
            vel_g[:, _] = interpolate.interpn(grid, vel_turb[_,:], pos_gas, method='linear')
            print(f'Dim {_+1} done.')

        # normalize so that RMS is 1
        vel_g = vel_g / np.sqrt(np.sum(vel_g**2, axis=1).mean())

        vel_g_bulge = vel_g * np.mean(sigma_turb) * 3**0.5
        vel_g_disc = vel_g * np.mean(sigma_turb_d) * 3**0.5

        vel_g_bulge_disc = vel_g_bulge * (rho/(rho+rho_d))[:, np.newaxis] + vel_g_disc * (rho_d/(rho+rho_d))[:, np.newaxis]

        vel_g_bh = vel_g * np.sqrt(G_cgs * smbh_mass / (r[:, np.newaxis]))/unit_velocity_cgs

        vel_gas += vel_g_bulge_disc + vel_g_bh


    if turb_on:
        spec_energy_g_bulge_disc = ism_temp * k_B  * (gamma-1)/ (mu * m_p) / unit_spec_energy_cgs
    else:
        spec_energy_g_bulge_disc = gas_spec_energy * (rho/(rho+rho_d)) + gas_spec_energy_d * (rho_d/(rho+rho_d))
        smbh_T = G_cgs*smbh_mass/r * mu*m_p/k_B
        gas_spec_energy_g_smbh = smbh_T * k_B  * (gamma-1)/ (mu * m_p) / unit_spec_energy_cgs
        spec_energy_g_bulge_disc += gas_spec_energy_g_smbh


    pos = pos_gas
    vel = vel_gas
    sp_egy_gas = np.zeros_like(pos_gas[:, 0]) + spec_energy_g_bulge_disc
    sp_egy = sp_egy_gas[:,np.newaxis]
    mass = mass_gas[:,np.newaxis]

    length = len(pos)
    pos = pos.astype(np.float64)
    vel = vel.astype(np.float64)
    sp_egy = sp_egy.astype(np.float32)
    mass = mass.astype(np.float32).ravel()

    del gas['Coordinates']
    gas.create_dataset('Coordinates', data=pos, compression="gzip")

    del gas['Velocities']
    gas.create_dataset('Velocities', data=vel, compression="gzip")

    gas.create_dataset('Masses', data=mass, compression="gzip")

    del gas['InternalEnergy']
    gas.create_dataset(
        'InternalEnergy', data=sp_egy.ravel(), compression="gzip")

    del gas['ParticleIDs']
    gas.create_dataset('ParticleIDs', data=(
        np.array(np.arange(length))+1).astype(np.uint), compression="gzip")


    del gas['Type']
    del gas['Density']

    if smbh_mass > 0:
        file.create_group('PartType5')
        file['PartType5'].create_dataset('Coordinates', data=[[0, 0, 0]])
        file['PartType5'].create_dataset(
            'ParticleIDs', data=[np.array(length+1).astype(np.uint)])
        file['PartType5'].create_dataset('Velocities', data=[[0, 0, 0]])
        file['PartType5'].create_dataset('Masses', data=[smbh_mass/ unit_mass_cgs], compression="gzip")

    header.attrs['NumPart_ThisFile'] = np.array((length, 0, 0, 0, 0, 0))
    header.attrs['NumPart_Total'] = np.array((length, 0, 0, 0, 0, 0))
    header.attrs['NumPart_Total_HighWord'] = np.array((0, 0, 0, 0, 0, 0))
    if smbh_mass > 0:
        header.attrs['NumPart_ThisFile'] = np.array((length, 0, 0, 0, 0, 1))
        header.attrs['NumPart_Total'] = np.array((length, 0, 0, 0, 0, 1))

    if turb_on:
        final_path = 'F://isothermal_ellipse_4M_turb_MT.hdf5'
    else:
        final_path = 'F://isothermal_ellipse_4M_noturb_MT.hdf5'

    # Reduce the final file size
    with h5py.File(final_path, 'w') as f2:
        for obj in file.keys():
            file.copy(obj, f2)
