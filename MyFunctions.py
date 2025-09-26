def TurbFieldML(res=257, kmin=2*np.pi, kmax=256*np.pi, nspec=-11/3, seed=1, method='cross', boxsize=1.0, ratio=1.0):
    '''
    In 3D Power spectrum P ~ k^(-11/3) is the equivalent to Kolmogorov energy spectrum E ~ k^(-5/3).
    Since P=A^2 -> A ~ k^(-13/6).
    For 'cross' method slope k^(-17/6) due to additional k for rotor in Fourier domain.
    Good read on ISM turbulence: arXiv:astro-ph/0404451v1 22 Apr 2004
    '''

    # Check for type
    method_list = ['cross', 'proj', 'f90']
    if method not in method_list:
        print('No such method:' + method)
        exit()
    if method == 'cross':
        slope = nspec-2
    if method == 'proj':
        slope = nspec
    if method == 'f90':
        slope = -nspec

    k = list(chain(range(0, -res//2, -1), range(res//2, 0, -1)))
    kvec = np.array(np.meshgrid(k, k, k))  # Generate cube res^3 in k space
    kmag = np.sqrt(np.sum(kvec**2, axis=0))
    kval = 2*np.pi/boxsize * kmag

    # Placeholder for potential
    Ak = np.zeros_like(kvec, dtype=complex)

    # Ak will be  0 outside ranges.
    mask = (kval >= kmin) & (kval <= kmax)

    rng = np.random
    state = rng.seed(seed)

    # Generates smooth power spectrum.
    pow_amp = 1.0

    # Get amplitude from Reileigh and phase from uniform and add variance to power spectrum

    if method == 'f90':
        kmin = 0.001
        powspec = pow_amp * (kval**2 + kmin**2)**(-(0.5*slope+1.))
        factor = np.sqrt(powspec * kval**(-0.5*slope+1.))

        ra1 = rng.uniform(seed, size=powspec.shape)
        ra2 = rng.uniform(seed, size=powspec.shape)
        gvar1 = np.cos(2 * np.pi * ra1) * np.sqrt(2)
        gvar2 = np.sin(2 * np.pi * ra1) * np.sqrt(2)

        amp_exp = np.sqrt(-np.log(ra2))

        gvar1 = gvar1 * amp_exp
        gvar2 = gvar2 * amp_exp

        for _ in range(3):  # For each dimmension
            Ak[_][mask] = factor * (gvar1 + 1j*gvar2)
    else:
        powspec = pow_amp * (kval[mask])**(slope)
        for _ in range(3):  # For each dimmension
            rand_reigh = rng.rayleigh(size=kval[mask].shape)
            rand_phase = rng.uniform(size=kval[mask].shape) * 2 * np.pi
            Ak[_][mask] = np.sqrt(powspec) * rand_reigh * np.exp(rand_phase*1j)

    '''
    At this point field is fully compressive.
    Converting compressional to solenoidal can be done in few ways.

    Rotor in Fourier space (cross):
    i*k x Ak - result should be divergence free.
    Real part of inverse fft is turbulent velocity field.
    Projection operator (proj):
    Alternative method is to use projection operator (arXiv:0808.0605; eq.1)
    A_i = (ratio*delta_ij + (1-2*ratio)*(k_i*k_j)/|k|^2)*A_j
    Final step is to do inverse fft of A_k and take only real part.
    '''

    if method == ('cross'):
        vel = fft.ifftn(np.cross(2.0*np.pi/boxsize*kvec*1j, Ak, axis=0)).real

    if method == 'proj':
        Ak_ = np.zeros_like(Ak)
        for i in range(3):
            for j in range(3):
                if i == j:
                    Ak_[i] += ratio * Ak[j]
                # np.divide handles division by 0 without warnings
                Ak_[i] += (1-2*ratio) * np.divide(kvec[i,]*kvec[j,], kval **
                                                  2, out=np.zeros_like(kval), where=kval != 0) * Ak[j]
        vel = fft.ifftn(Ak_).real

    return vel
