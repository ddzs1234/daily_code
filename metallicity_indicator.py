"""
Created on Sun Jan 26 22:16:49 2020

@author: ashley

D02-metallicity
KD02 
from Kewley 2008
"""
import sympy as sy
import numpy as np

def do2(xi, xi_e):
    """
    Parameters
    ----------
    xi : float array
        log10(NII_6583/Halpha)
    xi_e : float array
        error of xi.

    Returns
    -------
    float array
        metallicity, error

    """
    metal = []
    sigma = []
    for i in range(0, len(xi), 1):
        x = xi[i]
        x_e = xi_e[i]
        metal.append(9.12 + 0.73 * x)
        sigma.append(0.73 * x_e)
    return np.array(metal), np.array(sigma)



def kd02(logn2o, log23, logo32, logn2o_err, log23_err, logo32_err):
    """
    
    Parameters
    ----------
    logn2o : TYPE
        log10(NII/OII)
    log23 : TYPE
        log10((OII+OIII*4/3)/Hbeta)
    logo32 : TYPE
        log10(OIII/OII)
    logn2o_err : TYPE
        DESCRIPTION.
    log23_err : TYPE
        DESCRIPTION.
    logo32_err : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        metallicity and error

    """
    metal = []
    sigma = []
    for i in range(0, len(logn2o)):
        log_n2o = logn2o[i]
        log_23 = log23[i]
        log_o32 = logo32[i]
        n2o_e = logn2o_err[i]
        o32_e = logo32_err[i]
        r23_e = log23_err[i]

        if log_n2o > -1.2:
            z = sy.symbols('z', real=True)
            func = 1106.8660 - 532.15451 * z + 96.373260 * z**2 - \
                7.8106123 * z**3 + 0.23928247 * z**4 - log_n2o
            func1 = 1106.8660 - 532.15451 * z + 96.373260 * \
                z**2 - 7.8106123 * z**3 + 0.23928247 * z**4
            tmp = sy.solve(func, z)
            tmp_metal = [n for n in tmp if n > 8.4]
            if len(tmp_metal) >= 1:
                metal.append(np.min(tmp_metal))
                func1_dz = sy.diff(func1, z)
                dz = float(func1_dz.subs({z: np.min(tmp_metal)}))
                sigma_kd02 = n2o_e / dz
                sigma.append(sigma_kd02)

            elif len(tmp_metal) == 0:
                metal.append(np.max(tmp))
                func1_dz = sy.diff(func1, z)
                dz = float(func1_dz.subs({z: np.max(tmp)}))
                sigma_kd02 = n2o_e / dz
                sigma.append(sigma_kd02)
            else:
                print('wrong')
        else:
            # m91的金属丰度以及误差计算
            m91_a, m91_b = sy.symbols('m91_a m91_b')
            z_m91 = 12 - 4.944 + 0.767 * m91_a + 0.602 * m91_a**2 - m91_b * (
                0.29 + 0.332 * m91_a - 0.331 * m91_a**2)
            m91_da = sy.diff(z_m91, m91_a)
            m91_db = sy.diff(z_m91, m91_b)

            da = float(m91_da.subs({m91_a: log_23, m91_b: log_o32}))
            db = float(m91_db.subs({m91_a: log_23, m91_b: log_o32}))

            sigma_m91 = np.sqrt(da**2 * r23_e**2 + db**2 * o32_e**2)
            m91_z = float(z_m91.subs({m91_a: log_23, m91_b: log_o32}))

            # kk04金属丰度和误差计算
            kk04_a, kk04_b = sy.symbols('kk04_a kk04_b')

            q = ((32.81 - 1.153 * kk04_a**2 + 8.2 *
                  (-3.396 - 0.025 * kk04_a + 0.1444 * kk04_a**2)) /
                 (4.603 - 0.3119 * kk04_a - 0.163 * kk04_a**2 + 8.2 *
                  (-0.48 + 0.0271 * kk04_a + 0.02037 * kk04_a**2)))

            z_kk04 = 9.4 + 4.65 * kk04_b - 3.17 * kk04_b**2 - sy.log(q,10) * (
                0.272 + 0.547 * kk04_b - 0.513 * kk04_b**2)

            kk04_da = sy.diff(z_kk04, kk04_a)
            kk04_db = sy.diff(z_kk04, kk04_b)

            da1 = float(kk04_da.subs({kk04_a: log_o32, kk04_b: log_23}))
            db1 = float(kk04_db.subs({kk04_a: log_o32, kk04_b: log_23}))

            sigma_kk04 = np.sqrt(da1**2 * o32_e**2 + db1**2 * r23_e**2)
            kk04_z = float(z_kk04.subs({kk04_a: log_o32, kk04_b: log_23}))
            # kd02金属丰度
            sigma_kd02 = np.sqrt(sigma_m91**2 / 4 + sigma_kk04**2 / 4)
            sigma.append(sigma_kd02)

            metal.append((m91_z + kk04_z) / 2)
    return np.array(metal), np.array(sigma)   
