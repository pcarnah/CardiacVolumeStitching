import numpy as np

class MonogenicFilters:

    def __init__(self, xsize, ysize, zsize, wl_xy, wl_z):
        """
        Creates a set of frequency domain filters combining bandpass and Reisz
        filters for use in calculating the monogenic signal.

        :param xsize: X dimension
        :param yzise: Y dimension
        :param zsize: Z dimension
        :param wl_xy: Spatial wavelengths for x,y dims
        :param wl_z:  Spatial wavelenths for z dim
        """

        sigmaOnf = 0.5
        ratio = 0.98

        assert np.shape(wl_xy) == np.shape(wl_z)

        self.xsize = xsize
        self.ysize = ysize
        self.zsize = zsize

        self.numFilt = np.shape(wl_xy)[0]

        # Dimensions
        ymid = int(np.floor(self.ysize / 2))
        xmid = int(np.floor(self.xsize / 2))
        zmid = int(np.floor(self.zsize / 2))

        if (self.ysize % 2) == 0:
            ymax = ymid - 1
        else:
            ymax = ymid

        if (self.xsize % 2) == 0:
            xmax = xmid - 1
        else:
            xmax = xmid

        if (self.zsize % 2) == 0:
            zmax = zmid - 1
        else:
            zmax = zmid

        # Create the grid for the filter
        xGrid, yGrid, zGrid = np.meshgrid(range(-xmid, xmax + 1), range(-ymid, ymax + 1), range(-zmid, zmax + 1))
        xGrid = xGrid.transpose([1, 0, 2])
        yGrid = yGrid.transpose([1, 0, 2])
        zGrid = zGrid.transpose([1, 0, 2])

        xGrid = np.fft.ifftshift(xGrid)
        yGrid = np.fft.ifftshift(yGrid)
        zGrid = np.fft.ifftshift(zGrid)

        xGrid = xGrid / self.xsize
        yGrid = yGrid / self.ysize
        zGrid = zGrid / self.zsize

        self.bpFilt = np.zeros((self.xsize, self.ysize, self.zsize, self.numFilt))

        for i in range(self.numFilt):
            # Construct the filter -  first calculate the radial filter component
            f0 = 1.0 / wl_xy[i]  # Centre frequency of spatial filter.
            w0 = f0  # Normalised radius from centre of frequency plane corresponding to f0

            g0 = 1.0 / wl_z[i]  # Centre frequency of temporal filter.
            u0 = g0 / 0.5  # Normalised radius from centre of frequency plane corresponding to f0

            # Determine the spatial regions to use
            w = np.sqrt((xGrid ** 2) / (w0 ** 2) + (yGrid ** 2) / (w0 ** 2) + (zGrid ** 2) / (u0 ** 2))
            w[0, 0, 0] = 1  # Avoids division by zero

            # Computation of 3D log-gabor filter across the volumetric range
            self.bpFilt[:, :, :, i] = np.exp((-(np.log(w)) ** 2) / (2 * np.log(sigmaOnf) ** 2))

            # Set the DC value of the filter
            self.bpFilt[0, 0, 0, i] = 0

            # Also remove unwanted high frequency components in filters with even dimensions
            if self.xsize % 2 == 0:
                self.bpFilt[self.xsize // 2 + 1, :, 1, i] = 0
            if self.ysize % 2 == 0:
                self.bpFilt[:, self.ysize // 2 + 1, 1, i] = 0
            if self.zsize % 2 == 0:
                self.bpFilt[:, :, self.zsize // 2 + 1, i] = 0

        # Normalise by the maximum value of the sum of all filters
        sumFilt = np.sum(self.bpFilt, 3)
        self.bpFilt = self.bpFilt / np.max(sumFilt[:])

        # Generate the Riesz filter components (i.e. the odd filter whose components are imaginary)
        w = np.sqrt(yGrid ** 2 + xGrid ** 2 + zGrid ** 2)
        w[0, 0, 0] = 1
        self.ReiszFilt03 = 1 - (zGrid / w)
        self.ReiszFilt12 = (1j * xGrid - yGrid) / w

    def getMonogenicSignal(self, volume):
        return self.MonogenicSignal(volume, self)

    class MonogenicSignal:
        def __init__(self, volume, monogenic_filter):

            self.shape = (*volume.shape, monogenic_filter.numFilt)

            # Create output arrays
            self.Fm1 = np.zeros(self.shape)
            self.Fm2 = np.zeros(self.shape)
            self.Fm3 = np.zeros(self.shape)
            self.Fm4 = np.zeros(self.shape)

            # Compute the 3-dimensional fast Fourier transform of the original image
            F = np.fft.fftn(volume)

            # Filter using the Reisz filter
            R_03 = F * monogenic_filter.ReiszFilt03
            R_12 = F * monogenic_filter.ReiszFilt12

            # Compute the parts of the monogenic signal
            for flt in range(monogenic_filter.numFilt):
                F03 = np.fft.ifftn(R_03 * monogenic_filter.bpFilt[:, :, :, flt])
                F12 = np.fft.ifftn(R_12 * monogenic_filter.bpFilt[:, :, :, flt])
                self.Fm1[:, :, :, flt] = np.real(F03)
                self.Fm2[:, :, :, flt] = np.real(F12)
                self.Fm3[:, :, :, flt] = np.imag(F12)
                self.Fm4[:, :, :, flt] = np.imag(F03)

            self.odd = np.zeros(self.shape)
            np.sqrt(self.Fm2 ** 2 + self.Fm3 ** 2 + self.Fm4 ** 2, self.odd)

            self.even = np.abs(self.Fm1)

        def featureSymmetry(self, T=0.18):

            epsilon = np.finfo(float).eps

            odd = self.odd
            even = self.even

            # Calculate the denominator (= local energy + epsilon)
            denominator = np.sqrt(even ** 2 + odd ** 2) + epsilon

            # Calculate the numerators for FA and FS at all scales NB no need to take absolute of 'odd' as it must
            #  be positive due to the way it's calculated
            FS_numerator = np.clip(even - odd - T, 0, None)
            FA_numerator = np.clip(odd - even - T, 0, None)

            # Divide numerator by denominator
            FS = FS_numerator / denominator
            FA = FA_numerator / denominator

            # Sum across scales, and divide by number of scales to give value between 0
            # and 1 (i.e. take mean across the scale dimension)
            FS = np.mean(FS, 3)
            FA = np.mean(FA, 3)

            return FS, FA

        def localEnergy(self):
            LE = self.even ** 2 + self.odd ** 2
            return LE

        def localOrientation(self):
            epsilon = np.finfo(float).eps
            LO = np.zeros((3, *self.shape))

            denominator = self.odd + epsilon
            LO[0, :, :, :, :] = self.Fm2 / denominator
            LO[1, :, :, :, :] = self.Fm3 / denominator
            LO[2, :, :, :, :] = self.Fm4 / denominator
            return LO

        def localPhase(self):
            return np.arctan2(self.odd, self.Fm1)

        def orientedSymmetry(self, T=0.18):
            epsilon = np.finfo(float).eps

            # Combine the odd components
            odd = self.odd
            even = self.even

            # Calculate the denominator (= local energy + epsilon)
            denominator = np.sqrt(even ** 2 + odd ** 2) + epsilon

            # Calculate the numerators for FA and FS at all scales NB no need to take absolute of 'odd' as it must
            #  be positive due to the way it's calculated
            FS_numerator = np.clip(even - odd - T, 0, None)
            FA_numerator = np.clip(odd - even - T, 0, None)

            # Divide numerator by denominator and include orientation
            FA_x = (FA_numerator / denominator) * (self.Fm2 / (odd + epsilon))
            FA_y = (FA_numerator / denominator) * (self.Fm3 / (odd + epsilon))
            FA_z = (FA_numerator / denominator) * (self.Fm4 / (odd + epsilon))
            FS = (FS_numerator / denominator) * np.sign(self.Fm1)

            # Sum across scales, and divide by number of scales to give value between 0
            # and 1 (i.e. take mean across the scale dimension)
            FA_y = np.mean(FA_y, 3)
            FA_x = np.mean(FA_x, 3)
            FA_z = np.mean(FA_z, 3)
            FS = np.mean(FS, 3)

            return FA_x, FA_y, FA_z, FS


if __name__ == "__main__":
    shape = (224,144,208)
    filt = MonogenicFilters(*shape, [10, 20, 30, 40], [10, 20, 30, 40])
    mono = filt.getMonogenicSignal(np.zeros(shape))
    mono.featureSymmetry()
    mono.localEnergy()
    mono.localOrientation()
    mono.localPhase()
    _,_,_,FS = mono.orientedSymmetry()

    print(filt.bpFilt[0:10, 1, 1, 1])
    print(filt.ReiszFilt03[0:10, 1, 1])
    print(filt.ReiszFilt12[0:10, 1, 1])
