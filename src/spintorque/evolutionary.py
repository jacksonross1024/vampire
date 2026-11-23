from math import sqrt,sin,cos,exp 
import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

#   #------------------------------------------------
#   # Function to set transformation matrix m -> m'
#   #------------------------------------------------
def set_inverse_transformation_matrix(reference_vector, itm):

    # load vector into temporary variables for readability
    x = reference_vector[0]
    y = reference_vector[1]
    z = reference_vector[2]

    # Transform vectors in local coordinate system
    D1 = sqrt(y * y + z * z)
    D2 = sqrt(x * x + y * y + z * z)
    iD1 = 1.0 / D1
    iD2 = 1.0 / D2

    # define result matrix
    tm = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]

    tm[0][0] = D1 * iD2  # D1/D2
    tm[0][1] = -x * y * iD1 * iD2  # -x*y/(D1*D2)
    tm[0][2] = -x * z * iD1 * iD2  # -x*z/(D1*D2)

    tm[1][0] = 0.0
    tm[1][1] = z * iD1  # z/D1
    tm[1][2] = -y * iD1  # -y/D1;

    tm[2][0] = x * iD2  # x/D2
    tm[2][1] = y * iD2  # y/D2
    tm[2][2] = z * iD2  # z/D2

    # Set determinants for inverse matrix
    det1 = tm[0][0] * tm[1][1] * tm[2][2] + tm[0][1] * tm[1][2] * tm[2][0] + tm[0][2] * tm[1][0] * tm[2][1]
    det2 = tm[1][2] * tm[2][1] * tm[0][0] + tm[0][2] * tm[2][0] * tm[1][1] + tm[0][1] * tm[1][0] * tm[2][2]
    det_T = det1 - det2
    inv_det = 1.0 / det_T

    # Calculate inverse transformation matrix
    itm[0][0] = inv_det * (tm[1][1] * tm[2][2] - tm[1][2] * tm[2][1])
    itm[0][1] = inv_det * (tm[0][2] * tm[2][1] - tm[0][1] * tm[2][2])
    itm[0][2] = inv_det * (tm[0][1] * tm[1][2] - tm[0][2] * tm[1][1])
    itm[1][0] = inv_det * (tm[1][2] * tm[2][0] - tm[1][0] * tm[2][2])
    itm[1][1] = inv_det * (tm[0][0] * tm[2][2] - tm[0][2] * tm[2][0])
    itm[1][2] = inv_det * (tm[0][2] * tm[1][0] - tm[0][0] * tm[1][2])

    itm[2][0] = inv_det * (tm[1][0] * tm[2][1] - tm[1][1] * tm[2][0])
    itm[2][1] = inv_det * (tm[0][1] * tm[2][0] - tm[0][0] * tm[2][1])
    itm[2][2] = inv_det * (tm[0][0] * tm[1][1] - tm[0][1] * tm[1][0])
    itm[1][2] = inv_det * (tm[0][2] * tm[1][0] - tm[0][0] * tm[1][2])
    itm[2][0] = inv_det * (tm[1][0] * tm[2][1] - tm[1][1] * tm[2][0])
    itm[2][1] = inv_det * (tm[0][1] * tm[2][0] - tm[0][0] * tm[2][1])
    itm[2][2] = inv_det * (tm[0][0] * tm[1][1] - tm[0][1] * tm[1][0])

    return


#--------------------------------------------------------------------------------------------------
# Function to transform reference vector (rv) using transformation matrix (tm) moving rv -> result
#--------------------------------------------------------------------------------------------------
def transform_vector(rv, tm):

    # Declare result vector
    result = [0.0, 0.0, 0.0]

    # Calculate transformation
    result[0] = tm[0][0] * rv[0] + tm[0][1] * rv[1] + tm[0][2] * rv[2]
    result[1] = tm[1][0] * rv[0] + tm[1][1] * rv[1] + tm[1][2] * rv[2]
    result[2] = tm[2][0] * rv[0] + tm[2][1] * rv[1] + tm[2][2] * rv[2]

    return result


#------------------------------------------------------------------------------------------------
# Black magic code for Gaussian elimination
#------------------------------------------------------------------------------------------------

def gaussian_elimination(M, V):

    # initialise temporary arrays
    a_array = [
        [M[0][0], M[0][1], M[0][2]],
        [M[1][0], M[1][1], M[1][2]],
        [M[2][0], M[2][1], M[2][2]],
    ]
    b_array = [V[0], V[1], V[2]]
    # initialise result array
    x = [0.0, 0.0, 0.0]

    n = 3

    # Rearranging matrix(pivot)
    for i in range(n):  # column i
        tem1 = abs(a_array[i][i])
        p = i
        for j in range(i + 1, n):
            tem1 = abs(tem1)
            tem2 = abs(a_array[j][i])
            if tem2 > tem1:
                p = j
                tem1 = a_array[j][i]
        # row exchange in both the matrix
        for j in range(n):
            tem3 = a_array[i][j]
            a_array[i][j] = a_array[p][j]
            a_array[p][j] = tem3
        tem4 = b_array[i]
        b_array[i] = b_array[p]
        b_array[p] = tem4

    # Gaussian Elimination
    for i in range(n - 1):
        for j in range(i + 1, n):
            ratio = a_array[j][i] / a_array[i][i]  # diagonal value should not be zero
            for count in range(i, n):
                a_array[j][count] -= ratio * a_array[i][count]
            b_array[j] -= ratio * b_array[i]

    x[n - 1] = b_array[n - 1] / a_array[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        temp = b_array[i]
        for j in range(i + 1, n):
            temp -= a_array[i][j] * x[j]
        x[i] = temp / a_array[i][i]

    # copy to three vector
    result = [x[0], x[1], x[2]]

    return result


# material[2]:spin-diffusion-length = 36 #
# material[2]:spin-polarisation-conductivity = 0.23 #beta 
# material[2]:spin-polarisation-diffusion = 0.56 #
# material[2]:spin-accumulation = 1.48e7 #C/m^3 -> 0.11 e-/u.c.1.48e7 C/m^3 (Khmelevski) * 0.05 eV 
# material[2]:diffusion-constant = 0.002227 #m^2/s 
# material[2]:sd-exchange-constant = 0.05 !eV #50 meV 8.010883e-21 J, Oppeneer 

# material[2]:sot-spin-diffusion-length = 2.03 #
# material[2]:sot-spin-polarisation-conductivity = 0.31 #beta 
# material[2]:sot-spin-polarisation-diffusion = 0.31
# material[2]:sot-spin-accumulation = 1.48e7 
# material[2]:sot-diffusion-constant = 0.0007 
# material[2]:sot-sd-exchange-constant =  0.3 !eV 


# #Au values have minimal impact on Mn acc currently
# material[1]:spin-diffusion-length = 2000.0
# material[1]:spin-polarisation-conductivity = 0.99 
# material[1]:spin-polarisation-diffusion = 0.99     
# material[1]:spin-accumulation = 0.0# to match Mn acc when current on
# material[1]:diffusion-constant = 1e-3 #m^2/s 
# material[1]:sd-exchange-constant =  8.010883e-25 #5 ueV


# #Au values have minimal impact on Mn acc currently
# material[1]:sot-spin-diffusion-length = 2000.0
# material[1]:sot-spin-polarisation-conductivity = 0.99 
# material[1]:sot-spin-polarisation-diffusion = 0.99     
# material[1]:sot-spin-accumulation = 1.48# to match Mn acc when current on
# material[1]:sot-diffusion-constant = 0.0007 #m^2/s 
# material[1]:sot-sd-exchange-constant =  0.5 !eV

# def set_transfer_properties():
    # for cell in range(0, max_cells):

    #     sot_B  = sot_beta_cond[cell]
    #     sot_Bp = sot_beta_diff[cell]
    #     lambda_sdl = sot_lambda_sdl[cell]
    #     sot_Do = sot_diffusion[cell]
    #     sot_Jsd = sot_sd_exchange[cell]

    #     sot_BBp = 1.0/sqrt(1.0-sot_B*sot_Bp)
    #     sot_lambda_sf = lambda_sdl*sot_BBp
    #     sot_lambda_j = sqrt(2.0*hbar*sot_Do /sot_Jsd)
    #     sot_lambda_sf2 = sot_lambda_sf*sot_lambda_sf
    #     sot_lambda_j2 = sot_lambda_j*sot_lambda_j

    #     # l_perp = 0
    #     # l_L = 0
    #     sot_lambda_phi2 = 2*sot_lambda_j2
    #     sot_lambda_trans2 = sot_lambda_phi2 + sot_lambda_sf2
    #     #inside (1.0/lambda_trans2, -1.0/lambda_j2); 

    #     sot_inside = np.complex128(1.0/sot_lambda_sf2, -1.0/sot_lambda_trans2)
    #     sot_inv_lplus = np.sqrt([sot_inside])

    #     sot_a.append(np.real(sot_inv_lplus))
    #     sot_b.append(-np.imag(sot_inv_lplus))


def calculate_sot_accumulation(params, theta = np.pi*0.25, save=False):
    params = np.asarray(params, dtype=float)
    sot_a.clear()
    sot_b.clear()

    spin_galvanic_angle = params[0]
    #n_up(E_f)-n_dn(E_f)/(n_up(E_f)+n_dn(E_f))
    Au_sot_beta_cond = params[1]
    Mn_sot_beta_cond = params[2]
    sot_beta_cond = [Au_sot_beta_cond, Mn_sot_beta_cond, Mn_sot_beta_cond, Au_sot_beta_cond, Mn_sot_beta_cond, Mn_sot_beta_cond, Au_sot_beta_cond] ## beta

    #D_up(E_f)-D_dn(E_f)/(D_up(E_f)+D_dn(E_f))
    Au_sot_beta_diff = params[3]
    Mn_sot_beta_diff = params[4]
    sot_beta_diff = [Au_sot_beta_diff, Mn_sot_beta_diff, Mn_sot_beta_diff, Au_sot_beta_diff, Mn_sot_beta_diff, Mn_sot_beta_diff, Au_sot_beta_diff] ## beta_prime

    #s^2/m
    Au_sot_diff = params[5]
    Mn_sot_diff = params[6]
    sot_diffusion = [Au_sot_diff, Mn_sot_diff, Mn_sot_diff, Au_sot_diff, Mn_sot_diff, Mn_sot_diff, Au_sot_diff]

    sot_sa_source = [True, False, False, True, False, False, True] 

    #m
    Au_sot_lambda_sdl = params[7]*1e-10
    Mn_sot_lambda_sdl = params[8]*1e-10
    sot_lambda_sdl = [Au_sot_lambda_sdl, Mn_sot_lambda_sdl, Mn_sot_lambda_sdl, Au_sot_lambda_sdl, Mn_sot_lambda_sdl, Mn_sot_lambda_sdl, Au_sot_lambda_sdl]

    #C/m^3
    Au_sa = 1.48e7
    Mn_sa = 1.48e7
    sot_sa_inf = [Au_sa, Mn_sa, Mn_sa, Au_sa, Mn_sa, Mn_sa, Au_sa]

    #J
    Au_sot_sd_exchange = params[9]*1.60218e-19 
    Mn_sot_sd_exchange = params[10]*1.60218e-19
    sot_sd_exchange = [Au_sot_sd_exchange, Mn_sot_sd_exchange, Mn_sot_sd_exchange, Au_sot_sd_exchange, Mn_sot_sd_exchange, Mn_sot_sd_exchange, Au_sot_sd_exchange]


    m_init = np.zeros(7*3)

    for cell in range(0, max_cells):

        sot_B  = sot_beta_cond[cell]
        sot_Bp = sot_beta_diff[cell]
        lambda_sdl = sot_lambda_sdl[cell]
        sot_Do = sot_diffusion[cell]
        sot_Jsd = sot_sd_exchange[cell]

        sot_BBp = 1.0/sqrt(1.0-sot_B*sot_Bp)
        sot_lambda_sf = lambda_sdl*sot_BBp
        sot_lambda_j = sqrt(2.0*hbar*sot_Do /sot_Jsd)
        sot_lambda_sf2 = sot_lambda_sf*sot_lambda_sf
        sot_lambda_j2 = sot_lambda_j*sot_lambda_j

        # l_perp = 0
        # l_L = 0
        sot_lambda_phi2 = 2*sot_lambda_j2
        sot_lambda_trans2 = sot_lambda_phi2 + sot_lambda_sf2
        #inside (1.0/lambda_trans2, -1.0/lambda_j2); 

        sot_inside = np.complex128(1.0/sot_lambda_sf2, -1.0/sot_lambda_trans2)
        sot_inv_lplus = np.sqrt(sot_inside)

        sot_a.append(float(np.real(sot_inv_lplus)))
        sot_b.append(float(-np.imag(sot_inv_lplus)))

    for i in range(0,max_cells):
        mx = 0.0
        my = 0.0
        mz = 0.0
        if(i % 3 != 0):
            direction_sign = -1
            if i%3==2:
                direction_sign = 1 ##sign change for the up process
            
            mx = cos(theta)*direction_sign*sot_sa_inf[i]
            my = sin(theta)*direction_sign*sot_sa_inf[i]
            mz = 0.0

            m_init[3*i +0]  = cos(theta)*direction_sign
            m_init[3*i +1]  = sin(theta)*direction_sign
            m_init[3*i +2]  = mz

        sa_int[3*i + 0] = mx
        sa_int[3*i + 1] = my
        sa_int[3*i + 2] = mz

        
    itm = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0] ,[0.0, 0.0, 0.0]]  # inverse transformation matrix
    M = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0] , [0.0, 0.0, 0.0]]
    V = [0.0, 0.0, 0.0]  # general 3-vector

    # # reference basis vectors
    bx = [1.0,0.0,0.0]
    by = [0.0,1.0,0.0]
    bz = [0.0,0.0,1.0]

    # # local basis vectors
    b1 = [1.0,0.0,0.0]
    b2 = [0.0,1.0,0.0]
    b3 = [0.0,0.0,1.0]


    je_eff = fractional_electric_field_strength* je # current (C/s)

   
    i_muB = 1.0/9.274e-24 # # J/T
    i_e = 1.0/1.60217662e-19 # # electronic charge (Coulombs)
    microcell_volume = (micro_cell_size[0] *
                                    micro_cell_size[1] *
                                    micro_cell_thickness)*1.e-30 # # m^3
    atomcell_volume = 7.666015625e-30 # # muffin tin radii 15.7624e-30; #hard code for Mn2Au atom volume for now

    ##init
    for cell in range(0, max_cells):

        direction_sign = 1
        if (sot_sa_source[cell]):
            direction_sign = -1

        # # calculate cell id's
        cellx = 3*cell+0
        celly = 3*cell+1
        cellz = 3*cell+2

        pcell = cell - 1
        acell = cell + 1

        pcellx = 3*pcell + 0
        pcelly = 3*pcell + 1
        pcellz = 3*pcell + 2

        acellx = 3*acell + 0
        acelly = 3*acell + 1
        acellz = 3*acell + 2

        m_local = [0.0, 0.0, 0.0]
        # # copy array values to temporaries for readability
        if (sot_sa_source[cell]):
            m_local = [0.0,1.0,0.0] # Au cell artificial direction for correct basis of spin acc change
        else:
            m_local = (sa_int[cellx], sa_int[celly], sa_int[cellz]) # current cell magnetisations
        modm = sqrt(m_local[0]*m_local[0] + m_local[1]*m_local[1] + m_local[2]*m_local[2])
        # # Check for zero magnetization in normalization
        if modm > 1.e-11:
            m_local = (m_local[0]/modm, m_local[1]/modm, m_local[2]/modm)
        else:
            m_local = [0.0, 0.0, 0.0]
        # #---------------------------------------------------------------------
        # # Step 1 calculate coordinate transformation m -> m' for current cell
        # #---------------------------------------------------------------------

        # # Check for cell magnetization greater than 1e-8 mu_B
        if modm > 1.0e-11:
            # Initialise inverse transformation matrix
            set_inverse_transformation_matrix(m_local, itm)

            # Calculate basis vectors
            b1 = transform_vector(bz, itm)
            b2 = transform_vector(bx, itm)
            b3 = transform_vector(by, itm)
        else:
            # set default rotational frame
            b1 = bz
            b2 = bx
            b3 = by
        # #---------------------------------------------------------------------
        # # Step 2 determine coefficients for the spin accumulation
        # #---------------------------------------------------------------------

        # # Initialise temporary ants
        Bc = sot_beta_cond[cell] ## beta
        Bd = sot_beta_diff[cell] ## beta_prime
        Do = sot_diffusion[cell]
        jm0 = [0.0, 0.0, 0.0] ##different bounds for top and bottom layers
        if cell == 0:
            jm0 = (j_init_down_y[acellx],\
            j_init_down_y[acelly],\
            j_init_down_y[acellz])
        elif cell == max_cells-1:
            jm0 = (j_init_up_y[pcellx],\
            j_init_up_y[pcelly],\
            j_init_up_y[pcellz])
        else: 
            jm0 = (j_init_up_y[pcellx]+j_init_down_y[acellx],\
            j_init_up_y[pcelly]+j_init_down_y[acelly],\
            j_init_up_y[pcellz]+j_init_down_y[acellz])

        # #  Calculate gradient dsacc/dx
        # # need to include Roy's new interface dsa contribution for multilayers
        twoDo = 2.0*Do
        preD = twoDo*Bc*Bd

        M[0][0] = preD*m_local[0]*m_local[0] - twoDo     
        M[0][1] = preD*m_local[1]*m_local[0]             
        M[0][2] = preD*m_local[2]*m_local[0]
        M[1][0] = preD*m_local[0]*m_local[1]             
        M[1][1] = preD*m_local[1]*m_local[1] - twoDo     
        M[1][2] = preD*m_local[2]*m_local[1]
        M[2][0] = preD*m_local[0]*m_local[2]             
        M[2][1] = preD*m_local[1]*m_local[2]             
        M[2][2] = preD*m_local[2]*m_local[2] - twoDo

        # #no charge current conversion from out of plane process
        V[0] = jm0[0] ## - Bc*je_eff*m_local.x;
        V[1] = jm0[1] ## - Bc*je_eff*m_local.x;
        V[2] = jm0[2] ## - Bc*je_eff*m_local.z;

        divm_0 = gaussian_elimination(M, V)

        # # Calculate mp(0), c and d
        i_lsdl = 1.0/sot_lambda_sdl[cell]
        mp_inf = modm
        # #modify spin acc on Au depending on current. current Au m_inf = 0.0
        if(sot_sa_source[cell]):
            mp_inf *= fractional_electric_field_strength
        a_local = sot_a[cell]
        b_local = sot_b[cell]
        two_a = 2.0*a_local
        two_b = 2.0*b_local

        M[0][0] = -b1[0]*i_lsdl    
        M[0][1] = (-two_a*b2[0] + two_b*b3[0])    
        M[0][2] = (-two_b*b2[0] - two_a*b3[0])
        M[1][0] = -b1[1]*i_lsdl    
        M[1][1] = (-two_a*b2[1] + two_b*b3[1])    
        M[1][2] = (-two_b*b2[1] - two_a*b3[1])
        M[2][0] = -b1[2]*i_lsdl    
        M[2][1] = (-two_a*b2[2] + two_b*b3[2])    
        M[2][2] = (-two_b*b2[2] - two_a*b3[2])

        V[0] = divm_0[0] - b1[0]*mp_inf*i_lsdl
        V[1] = divm_0[1] - b1[1]*mp_inf*i_lsdl
        V[2] = divm_0[2] - b1[2]*mp_inf*i_lsdl
        coeff = gaussian_elimination(M, V)

        mp_0 = coeff[0]
        c    = coeff[1]
        d    = coeff[2]

        # #------------------------------------
        # # Step 3 calculate spin accumulation
        # #------------------------------------

        x = micro_cell_thickness*1.0e-10 ## Convert to metres
        cos_bx = cos(b_local*x)
        sin_bx = sin(b_local*x)
        e_xsdl = exp(-x*i_lsdl)
        e_ax   = exp(-a_local*x)
        prefac = (2.0*e_ax)

        sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl
        sa_perp2 = prefac*(c*cos_bx - d*sin_bx)
        sa_perp3 = prefac*(c*sin_bx + d*cos_bx)
                    
        # #--------------------------------------------
        # # Step 4 calculate the spin current (jm_end)
        # #--------------------------------------------
        ac_bd = a_local*c + b_local*d
        ad_bc = a_local*d - b_local*c

        divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl
        divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx)
        divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx)

        divsax = b1[0]*divsa_para + b2[0]*divsa_perp2 + b3[0]*divsa_perp3
        divsay = b1[1]*divsa_para + b2[1]*divsa_perp2 + b3[1]*divsa_perp3
        divsaz = b1[2]*divsa_para + b2[2]*divsa_perp2 + b3[2]*divsa_perp3

        dot = m_local[0]*divsax + m_local[1]*divsay + m_local[2]*divsaz
        pre_jmx = divsax - Bc*Bd*m_local[0]*dot
        pre_jmy = divsay - Bc*Bd*m_local[1]*dot
        pre_jmz = divsaz - Bc*Bd*m_local[2]*dot

        # #no charge current conversion
        jmx =   twoDo*pre_jmx
        jmy =   twoDo*pre_jmy
            # #modify spin current by current density and spin hall theta 
        if(sot_sa_source[cell]):
            jmy = je_eff*spin_galvanic_angle
        jmz =   twoDo*pre_jmz
        

        # #convert mp and m_perp
        sax = b1[0]*sa_para + b2[0]*sa_perp2 + b3[0]*sa_perp3
        say = b1[1]*sa_para + b2[1]*sa_perp2 + b3[1]*sa_perp3
        saz = b1[2]*sa_para + b2[2]*sa_perp2 + b3[2]*sa_perp3

        # # Save values. assume no 0 magnetisation cells
        sa_int[cellx] = sax
        sa_int[celly] = say
        sa_int[cellz] = saz

        # #init value saved in int for fast overwrite
        j_int_up_y[cellx] = jmx
        j_int_up_y[celly] = jmy
        j_int_up_y[cellz] = jmz

        # #down diffusion for Au spin current negative
        j_int_down_y[cellx] = jmx*direction_sign
        j_int_down_y[celly] = jmy*direction_sign
        j_int_down_y[cellz] = jmz*direction_sign
        
    ##up 1
    for cell in range(0, max_cells):
        direction_sign = 1 
        if cell%2==0:
            direction_sign = -1 ##sign change for the up process
        
        cellx = 3*cell+0
        celly = 3*cell+1
        cellz = 3*cell+2

        # # calculate previous cell id's
        pcellx = 3*(cell-1)+0
        pcelly = 3*(cell-1)+1
        pcellz = 3*(cell-1)+2

        # # calculate next cell id's
        acellx = 3*(cell+1)+0
        acelly = 3*(cell+1)+1
        acellz = 3*(cell+1)+2

        m_local = [0.0, 0.0, 0.0]
        # # copy array values to temporaries for readability
        if sot_sa_source[cell]:
            m_local[0] = 0
            m_local[1] = direction_sign # # Au cell artificial magnetisations
            m_local[2] = 0.0
        else:
            m_local[0] = sa_int[cellx]
            m_local[1] = sa_int[celly]
            m_local[2] = sa_int[cellz] # # current cell magnetisations

        modm = sqrt(m_local[0]*m_local[0] + m_local[1]*m_local[1] + m_local[2]*m_local[2])

        # # Check for zero magnetization in normalization
        if modm > 1.e-11:
            m_local[0] = m_local[0]/modm
            m_local[1] = m_local[1]/modm
            m_local[2] = m_local[2]/modm
        else:
            m_local[0] = 0.0
            m_local[1] = 0.0
            m_local[2] = 0.0

        #---------------------------------------------------------------------
        # Step 1 calculate coordinate transformation m -> m' for current cell
        #---------------------------------------------------------------------

        # Check for cell magnetization greater than 1e-8 mu_B
        if modm > 1.0e-11:
            # Initialise inverse transformation matrix
            set_inverse_transformation_matrix(m_local, itm)

            # Calculate basis vectors
            b1 = transform_vector(bz, itm)
            b2 = transform_vector(bx, itm)
            b3 = transform_vector(by, itm)
        else:
            # set default rotational frame
            b1 = bz
            b2 = bx
            b3 = by
        #---------------------------------------------------------------------
        # Step 2 determine coefficients for the spin accumulation
        #---------------------------------------------------------------------

        # Initialise temporary ants
        Bc = sot_beta_cond[cell] # beta
        Bd = sot_beta_diff[cell] # beta_prime
        Do = sot_diffusion[cell]
        jm0 = [0.0, 0.0, 0.0]
        #grab correct cell's init spin current (int label) bounds for top and bottom cells
        if direction_sign < 0 and cell != max_cells-1: 
            jm0[0] = j_int_down_y[acellx] 
            jm0[1] = j_int_down_y[acelly] 
            jm0[2] = j_int_down_y[acellz]
        elif cell != max_cells: 
            jm0[0] = j_int_up_y[pcellx] 
            jm0[1] = j_int_up_y[pcelly] 
            jm0[2] = j_int_up_y[pcellz]

        #  Calculate gradient dsacc/dx
        twoDo = 2.0*Do
        preD = twoDo*Bc*Bd

        M[0][0] = preD*m_local[0]*m_local[0] - twoDo     
        M[0][1] = preD*m_local[1]*m_local[0]             
        M[0][2] = preD*m_local[2]*m_local[0]
        M[1][0] = preD*m_local[0]*m_local[1]             
        M[1][1] = preD*m_local[1]*m_local[1] - twoDo     
        M[1][2] = preD*m_local[2]*m_local[1]
        M[2][0] = preD*m_local[0]*m_local[2]             
        M[2][1] = preD*m_local[1]*m_local[2]             
        M[2][2] = preD*m_local[2]*m_local[2] - twoDo

        #no charge current polarisation
        V[0] = jm0[0] # - Bc*je_eff*m_local.x;
        V[1] = jm0[1] # + Bc*je_eff*m_local.y;
        V[2] = jm0[2] # - Bc*je_eff*m_local.z;
        
        divm_0 = gaussian_elimination(M, V)

        # Calculate mp(0), c and d
        i_lsdl = 1.0/sot_lambda_sdl[cell]
        mp_inf = modm #sa_infinity[cell];
        if(sot_sa_source[cell]):
            mp_inf *= fractional_electric_field_strength
        a_local = sot_a[cell]
        b_local = sot_b[cell]
        two_a = 2.0*a_local
        two_b = 2.0*b_local

        M[0][0] = -b1[0]*i_lsdl    
        M[0][1] = (-two_a*b2[0] + two_b*b3[0])    
        M[0][2] = (-two_b*b2[0] - two_a*b3[0])
        M[1][0] = -b1[1]*i_lsdl    
        M[1][1] = (-two_a*b2[1] + two_b*b3[1])    
        M[1][2] = (-two_b*b2[1] - two_a*b3[1])
        M[2][0] = -b1[2]*i_lsdl    
        M[2][1] = (-two_a*b2[2] + two_b*b3[2])    
        M[2][2] = (-two_b*b2[2] - two_a*b3[2])

        V[0] = divm_0[0] - b1[0]*mp_inf*i_lsdl
        V[1] = divm_0[1] - b1[1]*mp_inf*i_lsdl
        V[2] = divm_0[2] - b1[2]*mp_inf*i_lsdl
        coeff = gaussian_elimination(M, V)

        mp_0 = coeff[0]
        c    = coeff[1]
        d    = coeff[2]

        #------------------------------------
        # Step 3 calculate spin accumulation
        #------------------------------------

        x = micro_cell_thickness*1.0e-10 # Convert to metres
        cos_bx = cos(b_local*x)
        sin_bx = sin(b_local*x)
        e_xsdl = exp(-x*i_lsdl)
        e_ax   = exp(-a_local*x)
        prefac = (2.0*e_ax)

        sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl
        sa_perp2 = prefac*(c*cos_bx - d*sin_bx)
        sa_perp3 = prefac*(c*sin_bx + d*cos_bx)
        
        #--------------------------------------------
        # Step 4 calculate the spin current (jm_end)
        #--------------------------------------------
        ac_bd = a_local*c + b_local*d
        ad_bc = a_local*d - b_local*c

        divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl
        divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx)
        divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx)

        divsax = b1[0]*divsa_para + b2[0]*divsa_perp2 + b3[0]*divsa_perp3
        divsay = b1[1]*divsa_para + b2[1]*divsa_perp2 + b3[1]*divsa_perp3
        divsaz = b1[2]*divsa_para + b2[2]*divsa_perp2 + b3[2]*divsa_perp3

        dot = m_local[0]*divsax + m_local[1]*divsay + m_local[2]*divsaz

        pre_jmx = divsax - Bc*Bd*m_local[0]*dot
        pre_jmy = divsay - Bc*Bd*m_local[1]*dot
        pre_jmz = divsaz - Bc*Bd*m_local[2]*dot
        jmx =   twoDo*pre_jmx
        jmy =   twoDo*pre_jmy
            #Au spin source sign direction dependent 
        if(sot_sa_source[cell]):
            jmy += direction_sign*je_eff*spin_galvanic_angle
        jmz =   twoDo*pre_jmz
        

        # Save values for the spin current
        if direction_sign < 0:
            j_int_down_y[cellx] = jmx
            j_int_down_y[celly] = jmy
            j_int_down_y[cellz] = jmz
        else:
            j_int_up_y[cellx] = jmx
            j_int_up_y[celly] = jmy
            j_int_up_y[cellz] = jmz 
    
    #int process down 1
    for cell in range(max_cells-1, -1, -1):
        direction_sign = -1 
        if cell%2==0:
            direction_sign = 1 #sign change for the down process

        cellx = 3*cell+0
        celly = 3*cell+1
        cellz = 3*cell+2

        # the sign change in the direction_sign will properly select the previous and next cells
        pcellx = 3*(cell-1)+0
        pcelly = 3*(cell-1)+1
        pcellz = 3*(cell-1)+2

        acellx = 3*(cell+1)+0
        acelly = 3*(cell+1)+1
        acellz = 3*(cell+1)+2

        m_local = [0.0, 0.0, 0.0]
        # copy array values to temporaries for readability
        if sot_sa_source[cell]:
            m_local = (0, direction_sign, 0.0) # Au cell artificial magnetisations
        else:
            m_local = (sa_int[cellx], sa_int[celly], sa_int[cellz]) # current cell magnetisations

        modm = sqrt(m_local[0]*m_local[0] + m_local[1]*m_local[1] + m_local[2]*m_local[2])

        # Check for zero magnetization in normalization
        if modm > 1.e-11:
            m_local = (m_local[0]/modm, m_local[1]/modm, m_local[2]/modm)
        else:
            m_local = [0.0, 0.0, 0.0]
        #--------------------------------------------------------------------
        # Step 1 calculate coordinate transformation m -> m' for current cell
        #---------------------------------------------------------------------

        # Check for cell magnetization greater than 1e-8 mu_B
        if modm > 1.0e-11:
            # Initialise inverse transformation matrix
            set_inverse_transformation_matrix(m_local, itm)

            # Calculate basis vectors
            b1 = transform_vector(bz, itm)
            b2 = transform_vector(bx, itm)
            b3 = transform_vector(by, itm)
        else:
            # set default rotational frame
            b1 = bz
            b2 = bx
            b3 = by
        #---------------------------------------------------------------------
        # Step 2 determine coefficients for the spin accumulation
        #---------------------------------------------------------------------

        # Initialise temporary ants
        Bc = sot_beta_cond[cell] # beta
        Bd = sot_beta_diff[cell] # beta_prime
        Do = sot_diffusion[cell]
        jm0 = [0.0, 0.0, 0.0]
        if direction_sign < 0 and cell != max_cells-1: 
            jm0 = (j_int_down_y[acellx], j_int_down_y[acelly], j_int_down_y[acellz])
        elif cell != 0:
            jm0 = (j_int_up_y[pcellx], j_int_up_y[pcelly], j_int_up_y[pcellz])

        #  Calculate gradient dsacc/dx
        twoDo = 2.0*Do
        preD = twoDo*Bc*Bd

        M[0][0] = preD*m_local[0]*m_local[0] - twoDo
        M[0][1] = preD*m_local[1]*m_local[0]
        M[0][2] = preD*m_local[2]*m_local[0]
        M[1][0] = preD*m_local[0]*m_local[1]
        M[1][1] = preD*m_local[1]*m_local[1] - twoDo
        M[1][2] = preD*m_local[2]*m_local[1]
        M[2][0] = preD*m_local[0]*m_local[2]
        M[2][1] = preD*m_local[1]*m_local[2]
        M[2][2] = preD*m_local[2]*m_local[2] - twoDo

        #no charge current conversion
        V[0] = jm0[0] # - Bc*je_eff*m_local.x;
        V[1] = jm0[1] # + Bc*je_eff*m_local.y;
        V[2] = jm0[2] # - Bc*je_eff*m_local.z;

        divm_0 = gaussian_elimination(M, V)

        # Calculate mp(0), c and d
        i_lsdl = 1.0/sot_lambda_sdl[cell]
        mp_inf = modm # sa_infinity[cell];
        if(sot_sa_source[cell]):
            mp_inf *= fractional_electric_field_strength
        a_local = sot_a[cell]
        b_local = sot_b[cell]
        two_a = 2.0*a_local
        two_b = 2.0*b_local

        M[0][0] = -b1[0]*i_lsdl
        M[0][1] = (-two_a*b2[0] + two_b*b3[0])
        M[0][2] = (-two_b*b2[0] - two_a*b3[0])
        M[1][0] = -b1[1]*i_lsdl
        M[1][1] = (-two_a*b2[1] + two_b*b3[1])
        M[1][2] = (-two_b*b2[1] - two_a*b3[1])
        M[1][2] = (-two_b*b2[1] - two_a*b3[1])
        M[2][0] = -b1[2]*i_lsdl
        M[2][1] = (-two_a*b2[2] + two_b*b3[2])
        M[2][2] = (-two_b*b2[2] - two_a*b3[2])

        V[0] = divm_0[0] - b1[0]*mp_inf*i_lsdl
        V[1] = divm_0[1] - b1[1]*mp_inf*i_lsdl
        V[2] = divm_0[2] - b1[2]*mp_inf*i_lsdl

        coeff = gaussian_elimination(M, V)

        mp_0 = coeff[0]
        c    = coeff[1]
        d    = coeff[2]

        #------------------------------------
        # Step 3 calculate spin accumulation
        #------------------------------------

        x = micro_cell_thickness*1.0e-10 # Convert to metres
        cos_bx = cos(b_local*x)
        sin_bx = sin(b_local*x)
        e_xsdl = exp(-x*i_lsdl)
        e_ax   = exp(-a_local*x)
        prefac = (2.0*e_ax)

        sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl
        sa_perp2 = prefac*(c*cos_bx - d*sin_bx)
        sa_perp3 = prefac*(c*sin_bx + d*cos_bx)
    
        #--------------------------------------------
        # Step 4 calculate the spin current (jm_end)
        #--------------------------------------------
        ac_bd = a_local*c + b_local*d
        ad_bc = a_local*d - b_local*c

        divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl
        divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx)
        divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx)

        divsax = b1[0]*divsa_para + b2[0]*divsa_perp2 + b3[0]*divsa_perp3
        divsay = b1[1]*divsa_para + b2[1]*divsa_perp2 + b3[1]*divsa_perp3
        divsaz = b1[2]*divsa_para + b2[2]*divsa_perp2 + b3[2]*divsa_perp3

        dot = m_local[0]*divsax + m_local[1]*divsay + m_local[2]*divsaz

        pre_jmx = divsax - Bc*Bd*m_local[0]*dot
        pre_jmy = divsay - Bc*Bd*m_local[1]*dot
        pre_jmz = divsaz - Bc*Bd*m_local[2]*dot
        jmx =   twoDo*pre_jmx
        jmy =   twoDo*pre_jmy
        if(sot_sa_source[cell]):
                jmy = direction_sign*je_eff*spin_galvanic_angle
        jmz =   twoDo*pre_jmz
        

        # Save values for the spin current
        if direction_sign < 0:
            j_int_down_y[cellx] = jmx
            j_int_down_y[celly] = jmy
            j_int_down_y[cellz] = jmz
        else:
            j_int_up_y[cellx] = jmx
            j_int_up_y[celly] = jmy
            j_int_up_y[cellz] = jmz

    #int process up, 2
    # second up sweep to use the newly calculated down sweep values
    for cell in range (0, max_cells):
        direction_sign = 1
        if cell%2==0 :
            direction_sign = -1
        
        cellx = 3*cell+0
        celly = 3*cell+1
        cellz = 3*cell+2

        # calculate previous cell id's
        pcellx = 3*(cell-1)+0
        pcelly = 3*(cell-1)+1
        pcellz = 3*(cell-1)+2

        acellx = 3*(cell+1)+0
        acelly = 3*(cell+1)+1
        acellz = 3*(cell+1)+2

        m_local = [0.0, 0.0, 0.0]
        # copy array values to temporaries for readability
        if sot_sa_source[cell]:
            m_local = (0, direction_sign, 0.0) # Au cell artificial magnetisations
        else:
            m_local = (sa_int[cellx], sa_int[celly], sa_int[cellz]) # current cell magnetisations

        modm = sqrt(m_local[0]*m_local[0] + m_local[1]*m_local[1] + m_local[2]*m_local[2])

        # Check for zero magnetization in normalization
        if modm > 1.e-11:
            m_local = (m_local[0]/modm, m_local[1]/modm, m_local[2]/modm)
        else:
            m_local = [0.0, 0.0, 0.0]
        #---------------------------------------------------------------------
        # Step 1 calculate coordinate transformation m -> m' for current cell
        #---------------------------------------------------------------------

        # Check for cell magnetization greater than 1e-8 mu_B
        if modm > 1.0e-11:
            # Initialise inverse transformation matrix
            set_inverse_transformation_matrix(m_local, itm)

            # Calculate basis vectors
            b1 = transform_vector(bz, itm)
            b2 = transform_vector(bx, itm)
            b3 = transform_vector(by, itm)
        else:
            # set default rotational frame
            b1 = bz
            b2 = bx
            b3 = by
        #---------------------------------------------------------------------
        # Step 2 determine coefficients for the spin accumulation
        #---------------------------------------------------------------------

        # Initialise temporary ants
        Bc = sot_beta_cond[cell] # beta
        Bd = sot_beta_diff[cell] # beta_prime
        Do = sot_diffusion[cell]
        jm0 = [0.0, 0.0, 0.0]
        if direction_sign < 0 and cell != max_cells-1: 
            jm0 = (j_int_down_y[acellx], j_int_down_y[acelly], j_int_down_y[acellz])
        elif cell != 0:
            jm0 = (j_int_up_y[pcellx], j_int_up_y[pcelly], j_int_up_y[pcellz])
        
        #  Calculate gradient dsacc/dx
        twoDo = 2.0*Do
        preD = twoDo*Bc*Bd

        M[0][0] = preD*m_local[0]*m_local[0] - twoDo
        M[0][1] = preD*m_local[1]*m_local[0]
        M[0][2] = preD*m_local[2]*m_local[0]
        M[1][0] = preD*m_local[0]*m_local[1]
        M[1][1] = preD*m_local[1]*m_local[1] - twoDo
        M[1][2] = preD*m_local[2]*m_local[1]
        M[2][0] = preD*m_local[0]*m_local[2]
        M[2][1] = preD*m_local[1]*m_local[2]
        M[2][2] = preD*m_local[2]*m_local[2] - twoDo

        # no charge current conversion
        V[0] = jm0[0] # - Bc*je_eff*m_local.x;
        V[1] = jm0[1] # + Bc*je_eff*m_local.y;
        V[2] = jm0[2] # - Bc*je_eff*m_local.z;

        divm_0 = gaussian_elimination(M, V)

        # Calculate mp(0), c and d
        i_lsdl = 1.0/sot_lambda_sdl[cell]
        mp_inf = modm # sa_infinity[cell];
        if(sot_sa_source[cell]):
            mp_inf *= fractional_electric_field_strength
        a_local = sot_a[cell]
        b_local = sot_b[cell]
        two_a = 2.0*a_local
        two_b = 2.0*b_local

        M[0][0] = -b1[0]*i_lsdl    
        M[0][1] = (-two_a*b2[0] + two_b*b3[0])    
        M[0][2] = (-two_b*b2[0] - two_a*b3[0])
        M[1][0] = -b1[1]*i_lsdl    
        M[1][1] = (-two_a*b2[1] + two_b*b3[1])    
        M[1][2] = (-two_b*b2[1] - two_a*b3[1])
        M[2][0] = -b1[2]*i_lsdl    
        M[2][1] = (-two_a*b2[2] + two_b*b3[2])    
        M[2][2] = (-two_b*b2[2] - two_a*b3[2])

        V[0] = divm_0[0] - b1[0]*mp_inf*i_lsdl
        V[1] = divm_0[1] - b1[1]*mp_inf*i_lsdl
        V[2] = divm_0[2] - b1[2]*mp_inf*i_lsdl
        coeff = gaussian_elimination(M, V)

        mp_0 = coeff[0]
        c    = coeff[1]
        d    = coeff[2]

        #------------------------------------
        # Step 3 calculate spin accumulation
        #------------------------------------

        x = micro_cell_thickness*1.0e-10 # Convert to metres
        cos_bx = cos(b_local*x)
        sin_bx = sin(b_local*x)
        e_xsdl = exp(-x*i_lsdl)
        e_ax   = exp(-a_local*x)
        prefac = (2.0*e_ax)

        sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl
        sa_perp2 = prefac*(c*cos_bx - d*sin_bx)
        sa_perp3 = prefac*(c*sin_bx + d*cos_bx)

        #--------------------------------------------
        # Step 4 calculate the spin current (jm_end)
        #--------------------------------------------
        ac_bd = a_local*c + b_local*d
        ad_bc = a_local*d - b_local*c

        divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl
        divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx)
        divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx)

        divsax = b1[0]*divsa_para + b2[0]*divsa_perp2 + b3[0]*divsa_perp3
        divsay = b1[1]*divsa_para + b2[1]*divsa_perp2 + b3[1]*divsa_perp3
        divsaz = b1[2]*divsa_para + b2[2]*divsa_perp2 + b3[2]*divsa_perp3

        dot = m_local[0]*divsax + m_local[1]*divsay + m_local[2]*divsaz

        pre_jmx = divsax - Bc*Bd*m_local[0]*dot
        pre_jmy = divsay - Bc*Bd*m_local[1]*dot
        pre_jmz = divsaz - Bc*Bd*m_local[2]*dot

        jmx =   twoDo*pre_jmx
        jmy =   twoDo*pre_jmy
        if(sot_sa_source[cell]):
            jmy += -direction_sign*je_eff*spin_galvanic_angle
        jmz =   twoDo*pre_jmz



        # Save values for the spin current
        if direction_sign < 0:
            j_int_down_y[cellx] = jmx
            j_int_down_y[celly] = jmy
            j_int_down_y[cellz] = jmz
        else:
            j_int_up_y[cellx] = jmx
            j_int_up_y[celly] = jmy
            j_int_up_y[cellz] = jmz
    
    
    #summation step
    for cell in range(0, max_cells):
        
        #direction sign for Au outgoing spin current for data output
        direction_sign = 1
        if sot_sa_source[cell]:
            direction_sign = -1
        cellx = 3*cell
        celly = 3*cell+1
        cellz = 3*cell+2

        m_local = [sa_int[cellx],  sa_int[celly],  sa_int[cellz]]
        modm = sqrt(m_local[0]*m_local[0] + m_local[1]*m_local[1] + m_local[2]*m_local[2])
        if modm > 1.e-11:
            m_local[0] = m_local[0]/modm
            m_local[1] = m_local[1]/modm
            m_local[2] = m_local[2]/modm
        else:
            m_local[0] = 0.0
            m_local[1] = 1.0 #default to +y for Au, get alternating sa in output but no physical change
            m_local[2] = 0.0
            modm = 1.0
        # Check for cell magnetization greater than 1e-8 mu_B
        if modm > 1.0e-11:
            # Initialise inverse transformation matrix
            set_inverse_transformation_matrix(m_local, itm)

            # Calculate basis vectors
            b1 = transform_vector(bz, itm)
            b2 = transform_vector(bx, itm)
            b3 = transform_vector(by, itm)
        else:
            # set default rotational frame
            b1 = bz
            b2 = bx
            b3 = by

        # Initialise temporary ants
        Bc = sot_beta_cond[cell] # beta
        Bd = sot_beta_diff[cell] # beta_prime
        Do = sot_diffusion[cell] # no contribution from x current since it was used to calculate the sa_int value at beginning of z axis steps
        jm0 = [j_int_up_y[cellx]+j_int_down_y[cellx], #+j_final_up_x[cellx], 
                            j_int_up_y[celly]+j_int_down_y[celly], #+j_final_up_x[cellx],
                            j_int_up_y[cellz]+j_int_down_y[cellz]] #+j_final_up_x[cellx]);                     

        
        twoDo = 2.0*Do
        preD = twoDo*Bc*Bd

        M[0][0] = preD*m_local[0]*m_local[0] - twoDo
        M[0][1] = preD*m_local[1]*m_local[0]
        M[0][2] = preD*m_local[2]*m_local[0]
        M[1][0] = preD*m_local[0]*m_local[1]
        M[1][1] = preD*m_local[1]*m_local[1] - twoDo
        M[1][2] = preD*m_local[2]*m_local[1]
        M[2][0] = preD*m_local[0]*m_local[2]
        M[2][1] = preD*m_local[1]*m_local[2]
        M[2][2] = preD*m_local[2]*m_local[2] - twoDo

        V[0] = jm0[0] # - Bc*je_eff*m_local.x;
        V[1] = jm0[1] #- Bc*je_eff*m_local.y;
        V[2] = jm0[2] # - Bc*je_eff*m_local.z;

        divm_0 = gaussian_elimination(M, V)

        i_lsdl = 1.0/sot_lambda_sdl[cell]
        mp_inf = modm # sa_infinity[cell];
        if(sot_sa_source[cell]):
            mp_inf *= fractional_electric_field_strength
        a_local = sot_a[cell]
        b_local = sot_b[cell]
        two_a = 2.0*a_local
        two_b = 2.0*b_local

        M[0][0] = -b1[0]*i_lsdl    
        M[0][1] = (-two_a*b2[0] + two_b*b3[0])    
        M[0][2] = (-two_b*b2[0] - two_a*b3[0])
        M[1][0] = -b1[1]*i_lsdl
        M[1][1] = (-two_a*b2[1] + two_b*b3[1])
        M[1][2] = (-two_b*b2[1] - two_a*b3[1])
        M[2][0] = -b1[2]*i_lsdl
        M[2][1] = (-two_a*b2[2] + two_b*b3[2])
        M[2][2] = (-two_b*b2[2] - two_a*b3[2])

        V[0] = divm_0[0] - b1[0]*mp_inf*i_lsdl
        V[1] = divm_0[1] - b1[1]*mp_inf*i_lsdl
        V[2] = divm_0[2] - b1[2]*mp_inf*i_lsdl

        coeff = gaussian_elimination(M, V)

        mp_0 = coeff[0]
        c    = coeff[1]
        d    = coeff[2]

        x = micro_cell_thickness*1.0e-10 # Convert to metres
        cos_bx = cos(b_local*x)
        sin_bx = sin(b_local*x)
        e_xsdl = exp(-x*i_lsdl)
        e_ax   = exp(-a_local*x)
        prefac = (2.0*e_ax)

        sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl
        sa_perp2 = prefac*(c*cos_bx - d*sin_bx)
        sa_perp3 = prefac*(c*sin_bx + d*cos_bx)

        #convert mp and m_perp
        sax = b1[0]*sa_para + b2[0]*sa_perp2 + b3[0]*sa_perp3
        say = b1[1]*sa_para + b2[1]*sa_perp2 + b3[1]*sa_perp3
        saz = b1[2]*sa_para + b2[2]*sa_perp2 + b3[2]*sa_perp3

        ac_bd = a_local*c + b_local*d
        ad_bc = a_local*d - b_local*c

        divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl
        divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx)
        divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx)

        divsax = b1[0]*divsa_para + b2[0]*divsa_perp2 + b3[0]*divsa_perp3
        divsay = b1[1]*divsa_para + b2[1]*divsa_perp2 + b3[1]*divsa_perp3
        divsaz = b1[2]*divsa_para + b2[2]*divsa_perp2 + b3[2]*divsa_perp3

        dot = m_local[0]*divsax + m_local[1]*divsay + m_local[2]*divsaz

        pre_jmx = divsax - Bc*Bd*m_local[0]*dot
        pre_jmy = divsay - Bc*Bd*m_local[1]*dot
        pre_jmz = divsaz - Bc*Bd*m_local[2]*dot

        jmx =   twoDo*pre_jmx
        jmy =   twoDo*pre_jmy
        if(sot_sa_source[cell]):
            jmy = je_eff*spin_galvanic_angle 
        jmz =   twoDo*pre_jmz

        sa_final[cellx] = sax
        sa_final[celly] = say
        sa_final[cellz] = saz
        

        prefac_sot = sot_sd_exchange[cell]*i_muB  #J T/J
        s = [m_init[3*cell], m_init[3*cell+1], m_init[3*cell+2]]
        mmod = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2])
        if(mmod > 0.0):
            s[0] /= mmod
            s[1] /= mmod
            s[2] /= mmod
        # spin_torque[cellx] = prefac_sot*(sax-s[0]*sot_sa_inf[cell])/sot_sa_inf[cell]
        # spin_torque[celly] = prefac_sot*(say-s[1]*sot_sa_inf[cell])/sot_sa_inf[cell]
        # spin_torque[cellz] = prefac_sot*(saz-s[2]*sot_sa_inf[cell])/sot_sa_inf[cell]
            # if(cell == 1) std::cout << "spin torque " <<  spin_torque[celly] << std::endl;
        spin_torque[cellx] = atomcell_volume * sot_sd_exchange[cell] * (sax) * i_e * i_muB #
        spin_torque[celly] = atomcell_volume * sot_sd_exchange[cell] * (say) * i_e * i_muB
        spin_torque[cellz] = atomcell_volume * sot_sd_exchange[cell] * (saz) * i_e * i_muB

        j_final_up_y[cellx] = jmx
        j_final_up_y[celly] = jmy
        j_final_up_y[cellz] = jmz

        j_final_down_y[cellx] = jmx*direction_sign
        j_final_down_y[celly] = jmy*direction_sign
        j_final_down_y[cellz] = jmz*direction_sign
    if(save):
        np.savetxt("sa-sot-solver-test.txt", np.stack( (sa_int, sa_final,  j_final_up_y, j_final_down_y, spin_torque ), axis=1), delimiter="\t", header="m, sa_final, d_sa,j_up, j_dn, ST")
    # return np.array([ #Sxs
    #     -(sa_int[7]*spin_torque[8] - sa_int[8]*spin_torque[7]),
    #     -(sa_int[8]*spin_torque[6] - sa_int[6]*spin_torque[8]),
    #     (sa_int[6]*spin_torque[7] - sa_int[7]*spin_torque[6]),
    # ])
    st = np.array([spin_torque[6], spin_torque[7], -spin_torque[8]])
    m_init = np.array([np.cos(theta), np.sin(theta), 0.0])*1e3
    torque = np.cross(m_init, st)
    return np.array([(sa_final[6]-sa_int[6])/sot_sa_inf[1], (sa_final[7]-sa_int[7])/sot_sa_inf[1], -(sa_final[8]-sa_int[8])/sot_sa_inf[1], torque[0], torque[1], torque[2]])

            

max_cells = 7
fractional_electric_field_strength = 1.0
je = 1e11 #A/m^2

micro_cell_size = (10.0, 10.0, 2.690333333)
micro_cell_thickness = 2.690333333

j_init_up_y = np.zeros(max_cells*3)
j_init_down_y = np.zeros(max_cells*3)
j_int_up_y = np.zeros(max_cells*3)
j_int_down_y = np.zeros(max_cells*3)

j_final_up_y = np.zeros(max_cells*3)
j_final_down_y = np.zeros(max_cells*3)

spin_torque = np.zeros(max_cells*3)

sa_final = np.zeros(max_cells*3)
sa_int = np.zeros(max_cells*3)

sot_a = []
sot_b = []
hbar = 1.05457162e-34

J_sd_ab = 1.0*1.60218e-19
mu_B = 9.274e-24

# set_transfer_properties()

theta_resolution = 72
ANGLES = np.linspace(0.01, np.pi-0.01, theta_resolution, endpoint=True)


def analytic_torque(theta, sd_ex, torque_scale):
    # torque_scale = 3.72e0
    # Example: replace with your actual analytic expression
    m = np.array([np.cos(theta), np.sin(theta), 0.0])
    s = [0.42*torque_scale*np.sin(2.0*theta),
          torque_scale*0.795*(np.sin(theta)*np.sin(theta) + 1.32),
         -0.12*torque_scale*np.cos(theta)]
    sa = [0.42*torque_scale*np.sin(2.0*theta)*sd_ex*1.60218e-19/mu_B,
          torque_scale*0.795*(np.sin(theta)*np.sin(theta) + 1.32)*sd_ex*1.60218e-19/mu_B,
         -0.12*torque_scale*np.cos(theta)*sd_ex*1.60218e-19/mu_B]
    
    # cross = np.cross(m,sa)
    cross = np.array([s[0], s[1], s[2], m[1]*sa[2]-m[2]*sa[1], m[2]*sa[0]-m[0]*sa[2], m[0]*sa[1] - m[1]*sa[0]])
    return cross


def simulate_torque(params, theta):
    """
    params: dict or array of optimization variables
    theta: angle on unit circle
    returns: torque vector (Tx, Ty, Tz)
    """
    # Your full simulation here
    # e.g. micromagnetics, FEM, numerical solver, etc.
    # NO gradients, NO assumptions
    torque = calculate_sot_accumulation(params, theta)
    return torque

def normalize(v):
    return v / (np.linalg.norm(v) + 1e-12)

def objective_function(params):
    params = np.asarray(params, dtype=float)
    errors = []

    for theta in ANGLES:
        T_sim = simulate_torque(params, theta)
        T_ref = analytic_torque(theta, params[10], 0.32*0.01)
       

        T_sim_l = T_sim[0]*T_sim[0]+T_sim[1]*T_sim[1]+T_sim[2]*T_sim[2]
        T_ref_l = T_ref[0]*T_ref[0]+T_ref[1]*T_ref[1]+T_ref[2]*T_ref[2] 

        T_sim_t_l = T_sim[3]*T_sim[3]+T_sim[4]*T_sim[4]+T_sim[5]*T_sim[5]
        T_ref_t_l = T_ref[3]*T_ref[3]+T_ref[4]*T_ref[4]+T_ref[5]*T_ref[5]
        loss =  (T_sim_l-T_ref_l)**2\
              + (T_sim_t_l-T_ref_t_l)**2\
              + 2*(1 - np.dot(normalize(T_sim[:2]), normalize(T_ref[:2])))\
              + sqrt(2/10.)*(T_sim[0] - T_ref[0])**2\
              + sqrt(2/10.)*(T_sim[1] - T_ref[1])**2\
              + sqrt(1.0/10.)*(T_sim[2] - T_ref[2])**2\
              + sqrt(2/10.)*(T_sim[3] - T_ref[3])**2\
              + sqrt(2/10.)*(T_sim[4] - T_ref[4])**2\
              + sqrt(1.0/10.)*(T_sim[5] - T_ref[5])**2
        
        # errors.append(loss)
        # err = 1 - np.dot(normalize(T_sim), normalize(T_ref))
        # err = 1 - np.dot(T_sim, T_ref)
        errors.append(loss)

    return np.mean(errors)


def init_params(starting_params = None):
    
    if(starting_params == None):
        starting_params = [2.83556557e-01, 9.81771826e-01, 3.45758999e-01, 8.97133929e-01,
        9.87768730e-01, 1.84952513e-05, 2.24134454e-04, 1.93123177e+03,
        2.89335349e+01, 2.58197851e-01, 1.77273585e-01]

    parametrization = ng.p.Array(init=starting_params).set_bounds(
        #      SGA  Au_Beta  Mn_B Au_B'   Mn_B' Au_D  Mn_D Au_lsf Mn_lsf Au_Jsd Mn_Jsd, torque_scaling
        lower=[0.15, 0.001, 0.25, 0.001, 0.001, 1e-5, 1e-5, 500.0,  2.5, 0.1, 0.1],
        upper=[0.29,  0.99,  0.35, 0.99,  0.99,  1e-3, 1e-3, 3000.0, 30.0, 1.2, 1.2]
    )
    return parametrization



optimize = True
seed_count = 5
seed = 313468

N_WORKERS = 8

def evaluate(candidate):
    return objective_function(candidate.value)


from concurrent.futures import ProcessPoolExecutor, as_completed

import nevergrad as ng 
if(optimize):

    N_BATCH = 11*11*20
    best_params = None
    best_loss = np.inf
    for seeded in range(seed_count):
        print("Running seed ", seeded)
        np.random.seed(seed + seeded)
        
        parametrization = init_params()
        parametrization.random_state.seed(seed+seeded)
        scale = 11*11*50

        optimizer = ng.optimizers.DE(parametrization=parametrization, budget=scale, num_workers=N_WORKERS)

        with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
            remaining = optimizer.budget
            while remaining > 0:
                batch = min(N_BATCH, remaining)
                candidates = [optimizer.ask() for _ in range(batch)]
                futures = {executor.submit(evaluate, c): c for c in candidates}

                for f in as_completed(futures):
                    c = futures[f]
                    optimizer.tell(c, f.result())
                    remaining -= 1

        recommendation = optimizer.provide_recommendation()

        if recommendation.loss < best_loss:
            best_params = recommendation.value
            best_loss = recommendation.loss


best_params = [2.53068891e-01, 9.12106757e-01, 3.12152876e-01, 6.40125584e-01,
9.13164158e-01, 8.16509160e-04, 5.00845174e-04, 2.20427293e+03,
1.66846306e+01, 2.17847606e-01, 1.01295724e-01]
parametrization = init_params(best_params)
print(parametrization.value)

scale = 11*11*50
N_BATCH = 11*11
optimizer = ng.optimizers.CMA(parametrization=parametrization, budget=scale,num_workers=N_WORKERS)

with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        remaining = optimizer.budget
        while remaining > 0:
            batch = min(N_BATCH, remaining)
            candidates = [optimizer.ask() for _ in range(batch)]
            futures = {executor.submit(evaluate, c): c for c in candidates}

            for f in as_completed(futures):
                c = futures[f]
                optimizer.tell(c, f.result())
                remaining -= 1

recommendation = optimizer.provide_recommendation()
best_params = recommendation.value
best_loss = recommendation.loss

print(best_params)



if plt is not None:
    angles = np.linspace(0.01, np.pi-0.01, 18, endpoint=True)

    T_sim = np.array([simulate_torque(best_params, th) for th in angles])
    T_ref = np.array([analytic_torque(th, best_params[10], 0.32*0.01) for th in angles])

    labels = ['dS_x/S', 'dS_y/S', 'dS_z/S', 'T_x', 'T_y', 'T_z']

    fig, axs = plt.subplots(2, 3, sharex=True, figsize=(8, 6))
    axs = np.atleast_1d(axs).ravel()

    for i in range(6):
        axs[i].plot(angles, (T_sim[:, i]), 'r', label='Simulated')
        axs[i].plot(angles, T_ref[:, i], 'k--', label='Analytic')
        
        axs[i].set_ylabel(labels[i])
        axs[i].legend()

    axs[-1].set_xlabel(r'$\theta$')
    plt.tight_layout()
    plt.show()

    mag_sim = np.linalg.norm(T_sim, axis=1)
    mag_ref = np.linalg.norm(T_ref, axis=1)

    plt.plot(angles, mag_ref, 'k--', label='|T_ref|')
    plt.plot(angles, mag_sim, 'r', label='|T_sim|')
    plt.xlabel(r'$\theta$')
    plt.ylabel('Torque magnitude')
    plt.legend()
    plt.show()

    angle_error = [
        1 - np.dot(normalize(T_sim[i]), normalize(T_ref[i]))
        for i in range(len(angles))
    ]

    plt.plot(angles, angle_error)
    plt.xlabel(r'$\theta$')
    plt.ylabel('1 - cos(angle)')
    plt.title('Directional error')
    plt.show()


calculate_sot_accumulation(best_params, save=True  )
