#  /*====================================================================
#  Copyright (C) 2019  Patrick Carnahan <pcarnah@uwo.ca>
#  Unauthorized copying of this file, via any medium is strictly prohibited
#  Proprietary and confidential
#  ====================================================================*/

import numpy as np
import numpy.fft as fft
import torch
from timeit import default_timer as timer

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

class MonogenicFilters:

    def __init__(self, size, spacing, wl_xy):
        """
        Creates a set of frequency domain filters combining bandpass and Reisz
        filters for use in calculating the monogenic signal.

        :param xsize: X dimension
        :param yzise: Y dimension
        :param zsize: Z dimension
        :param wl_xy: Spatial wavelengths for x,y dims
        """

        sigmaOnf = 0.5
        ratio = 0.98

        self.xsize = size[0]
        self.ysize = size[1]
        self.zsize = size[2]

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

        XR = torch.empty(self.xsize, dtype=torch.long, device=device)
        XR[:xmax + 1] = torch.arange(xmax + 1)
        XR[xmax + 1:] = torch.arange(-xmid, 0)

        YR = torch.empty(self.ysize, dtype=torch.long, device=device)
        YR[:ymax + 1] = torch.arange(ymax + 1)
        YR[ymax + 1:] = torch.arange(-ymid, 0)

        ZR = torch.empty(self.zsize, dtype=torch.long, device=device)
        ZR[:zmax + 1] = torch.arange(zmax + 1)
        ZR[zmax + 1:] = torch.arange(-zmid, 0)


        xGrid, yGrid, zGrid = torch.meshgrid([XR, YR, ZR])

        xGrid = xGrid.type(torch.float64) / self.xsize
        yGrid = yGrid.type(torch.float64) / self.ysize
        zGrid = zGrid.type(torch.float64) / self.zsize

        xGrid2 = (xGrid ** 2)
        yGrid2 = (yGrid ** 2)
        zGrid2 = (zGrid ** 2)

        # Construct the filter -  first calculate the radial filter component
        w0 = 1.0 / torch.as_tensor(wl_xy, dtype=torch.float64, device=device)  # Centre frequency of spatial filter.


        # Determine the spatial regions to use
        w = torch.sqrt(xGrid2.unsqueeze(len(size)).expand((*xGrid.shape, self.numFilt)) / ((w0 / spacing[0]) ** 2)
                       + yGrid2.unsqueeze(len(size)).expand((*xGrid.shape, self.numFilt)) / ((w0 / spacing[1]) ** 2)
                       + zGrid2.unsqueeze(len(size)).expand((*xGrid.shape, self.numFilt)) / ((w0 / spacing[2]) ** 2))
        w[0, 0, 0, :] = 1  # Avoids division by zero

        # Computation of 3D log-gabor filter across the volumetric range
        self.bpFilt = torch.exp((-(torch.log(w)) ** 2) / (2 * np.log(sigmaOnf) ** 2))


        # Set the DC value of the filter
        self.bpFilt[0, 0, 0, :] = 0

        # Also remove unwanted high frequency components in filters with even dimensions
        if self.xsize % 2 == 0:
            self.bpFilt[self.xsize // 2 + 1, :, 1, :] = 0
        if self.ysize % 2 == 0:
            self.bpFilt[:, self.ysize // 2 + 1, 1, :] = 0

        # Normalise by the maximum value of the sum of all filters
        sumFilt = torch.sum(self.bpFilt, 3)
        self.bpFilt = self.bpFilt / torch.max(sumFilt[:])
        self.bpFilt = torch.nn.functional.pad(self.bpFilt[:, :, :, :, None], [0, 1], "constant", 0)


        # Generate the Riesz filter components (i.e. the odd filter whose components are imaginary)
        w = torch.sqrt(xGrid2 + yGrid2 + zGrid2)
        w[0, 0, 0] = 1

        self.ReiszFilt03 = torch.zeros((*xGrid.shape, 2), dtype=torch.float64, device=device)
        self.ReiszFilt12 = torch.zeros((*xGrid.shape, 2), dtype=torch.float64, device=device)

        self.ReiszFilt03[:,:,:,0] = 1 - (zGrid / w)
        self.ReiszFilt12[:, :, :, 1] = xGrid / w
        self.ReiszFilt12[:, :, :, 0] = -yGrid / w


    def getMonogenicSignal(self, volume):
        return self.MonogenicSignal(volume, self)

    class MonogenicSignal:
        def __init__(self, volume, monogenic_filter):

            self.shape = (*volume.shape, monogenic_filter.numFilt)

            # Create output arrays
            self.Fm1 = torch.zeros(self.shape, dtype=torch.float32, device=device)
            self.Fm2 = torch.zeros(self.shape, dtype=torch.float32, device=device)
            self.Fm3 = torch.zeros(self.shape, dtype=torch.float32, device=device)
            self.Fm4 = torch.zeros(self.shape, dtype=torch.float32, device=device)

            # Compute the 3-dimensional fast Fourier transform of the original image
            V = torch.as_tensor(volume, device=device)
            VC = torch.nn.functional.pad(V.unsqueeze(3), [0, 1], "constant", 0)
            F = torch.fft(VC, 3)

            # Filter using the Reisz filter
            def complex_mult(A,B):
                R = torch.empty_like(A)

                R[:,:,:,0] = A.select(-1,0) * B.select(-1,0) - A.select(-1,1) * B.select(-1,1)
                R[:,:,:,1] = A.select(-1,1) * B.select(-1,0) + A.select(-1,0) * B.select(-1,1)
                return R


            R_03 = complex_mult(F, monogenic_filter.ReiszFilt03)
            R_12 = complex_mult(F, monogenic_filter.ReiszFilt12)

            # Compute the parts of the monogenic signal
            for flt in range(monogenic_filter.numFilt):
                t1 = complex_mult(R_03, monogenic_filter.bpFilt[:, :, :, flt, :])
                t2 = complex_mult(R_12, monogenic_filter.bpFilt[:, :, :, flt, :])
                F03 = torch.ifft(t1, 3)
                F12 = torch.ifft(t2, 3)
                self.Fm1[:, :, :, flt] = F03[:,:,:,0]
                self.Fm2[:, :, :, flt] = F12[:,:,:,0]
                self.Fm3[:, :, :, flt] = F12[:,:,:,1]
                self.Fm4[:, :, :, flt] = F03[:,:,:,1]

            self.odd = torch.sqrt(self.Fm2 ** 2 + self.Fm3 ** 2 + self.Fm4 ** 2)

            self.even = torch.abs(self.Fm1)

        def featureSymmetry(self, T=0.18):

            epsilon = np.finfo(float).eps

            odd = self.odd
            even = self.even

            # Calculate the denominator (= local energy + epsilon)
            denominator = torch.sqrt(even ** 2 + odd ** 2) + epsilon

            # Calculate the numerators for FA and FS at all scales NB no need to take absolute of 'odd' as it must
            #  be positive due to the way it's calculated
            FS_numerator = torch.clamp(even - odd - T, 0, float('inf'))
            FA_numerator = torch.clamp(odd - even - T, 0, float('inf'))

            # Divide numerator by denominator
            FS = FS_numerator / denominator
            FA = FA_numerator / denominator

            # Sum across scales, and divide by number of scales to give value between 0
            # and 1 (i.e. take mean across the scale dimension)
            FS = torch.mean(FS, 3)
            FA = torch.mean(FA, 3)

            return np.array(FS.cpu()), np.array(FA.cpu())

        def localEnergy(self):
            LE = self.even ** 2 + self.odd ** 2
            return np.array(LE.cpu())

        def localOrientation(self):
            epsilon = np.finfo(float).eps
            LO = torch.empty((3, *self.shape), dtype=torch.float32, device=device)

            denominator = self.odd + epsilon
            LO[0, :, :, :, :] = self.Fm2 / denominator
            LO[1, :, :, :, :] = self.Fm3 / denominator
            LO[2, :, :, :, :] = self.Fm4 / denominator
            return np.array(LO.cpu())

        def localPhase(self):
            return np.array(torch.atan2(self.odd, self.Fm1).cpu())

        def orientedSymmetry(self, T=0.18):
            epsilon = np.finfo(float).eps

            # Combine the odd components
            odd = self.odd
            even = self.even

            # Calculate the denominator (= local energy + epsilon)
            denominator = torch.sqrt(even ** 2 + odd ** 2) + epsilon

            # Calculate the numerators for FA and FS at all scales NB no need to take absolute of 'odd' as it must
            #  be positive due to the way it's calculated
            FS_numerator = torch.clamp(even - odd - T, 0, float('inf'))
            FA_numerator = torch.clamp(odd - even - T, 0, float('inf'))

            # Divide numerator by denominator and include orientation
            FA_x = (FA_numerator / denominator) * (self.Fm2 / (odd + epsilon))
            FA_y = (FA_numerator / denominator) * (self.Fm3 / (odd + epsilon))
            FA_z = (FA_numerator / denominator) * (self.Fm4 / (odd + epsilon))
            FS = (FS_numerator / denominator) * torch.sign(self.Fm1)

            # Sum across scales, and divide by number of scales to give value between 0
            # and 1 (i.e. take mean across the scale dimension)
            FA_y = torch.mean(FA_y, 3)
            FA_x = torch.mean(FA_x, 3)
            FA_z = torch.mean(FA_z, 3)
            FS = torch.mean(FS, 3)

            return np.array(FA_x.cpu()), np.array(FA_y.cpu()), np.array(FA_z.cpu()), np.array(FS.cpu())


if __name__ == "__main__":
    np.random.seed(10)
    a = np.random.rand(200,200,200)

    shape = a.shape
    spacing = [1, 1, 1]
    start = timer()
    filt = MonogenicFilters(shape, spacing, [10, 20, 30, 40, 50])
    end = timer()
    print("Time: {}".format(end - start))

    start = timer()
    mono = filt.getMonogenicSignal(a)
    end = timer()
    print("Time: {}".format(end - start))

    start = timer()
    mono.featureSymmetry()
    mono.localEnergy()
    mono.localOrientation()
    mono.localPhase()
    _,_,_,FS = mono.orientedSymmetry()

    end = timer()
    print("Time: {}".format(end - start))

    print(mono.Fm1[0:10, 1, 1, 1])
    print(mono.Fm2[0:10, 1, 1])
    print(mono.Fm3[0:10, 1, 1])
    print(mono.Fm4[0:10, 1, 1])

