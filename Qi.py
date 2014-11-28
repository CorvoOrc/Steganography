__author__ = 'Steshenko'

import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
import matplotlib.pyplot as plt
import matplotlib
import pywt, numpy
import pylab
import math

def embeding_subband(data, w_b, w_e, h_b, h_e, coeff):
    k, l = 0, 0
    for i in xrange(h_b, h_e):
        for j in xrange(w_b, w_e):
            data[i][j], l = coeff[k][l], l + 1
        k, l = k + 1, 0
        
    return data

def draw_edge(data, w, h, level):
    W, H = w, h
    for k in xrange(4):
        for i in xrange(W):
            data[H/2][i] = -1
        H, W = H/2, W/2
        
    W, H = w, h
    for k in xrange(4):
        for i in xrange(H):
            data[i][W/2] = -1
        H, W = H/2, W/2

    return data

def embeding_subbands(w, h, coeffs, level):
    LL4, (LH4, HL4, HH4), (LH3, HL3, HH3), (LH2, HL2, HH2), (LH1, HL1, HH1) = coeffs
    matrixCoeffs = [[0 for j in xrange(w)] for i in xrange(h)]
    
    matrixCoeffs = embeding_subband(matrixCoeffs, 0,    w/16, 0,    h/16, LL4)
    matrixCoeffs = embeding_subband(matrixCoeffs, 0,    w/16, h/16, h/8,  LH4)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/16, w/8,  0,    h/16, HL4)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/16, w/8,  h/16, h/8,  HH4)
    matrixCoeffs = embeding_subband(matrixCoeffs, 0,    w/8,  h/8,  h/4,  LH3)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/8,  w/4,  0,    h/8,  HL3)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/8,  w/4,  h/8,  h/4,  HH3)
    matrixCoeffs = embeding_subband(matrixCoeffs, 0,    w/4,  h/4,  h/2,  LH2)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/4,  w/2,  0,    h/4,  HL2)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/4,  w/2,  h/4,  h/2,  HH2)
    matrixCoeffs = embeding_subband(matrixCoeffs, 0,    w/2,  h/2,  h,    LH1)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/2,  w,    0,    h/2,  HL1)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/2,  w,    h/2,  h,    HH1)

    matrixCoeffs = draw_edge(matrixCoeffs, w, h, level)
    
    return matrixCoeffs
    
	
def marking_rand_matr(data, w_wm, h_wm):
    random_matrix = [[0 for j in xrange(w_wm)] for i in xrange(h_wm)]
    for i in xrange(h_wm):
        for j in xrange(w_wm):
            if data.getpixel((j, i)) == 255: random_matrix[i][j] = 1
            elif data.getpixel((j, i)) == 0: random_matrix[i][j] = -1

    return random_matrix

def varience(x1, x2, x3, x4, n):
    return 1/n*(x1**2 + x2**2 + x3**2 + x4**2) - (1/n*(x1 + x2 + x3 + x4))**2
    
def watermarked(OriginalSubbandValue, w, TotalChange, PseudorandomSequence):
    WatermarkedSubbandValue = TotalChange
    
    for i in xrange(len(PseudorandomSequence)):
        for j in xrange(len(PseudorandomSequence[i])):
            WatermarkedSubbandValue[i][j] *= PseudorandomSequence[i][j]
    
    for i in xrange(len(WatermarkedSubbandValue)):
        for j in xrange(len(WatermarkedSubbandValue[i])):
            WatermarkedSubbandValue[i][j] *= w

    for i in xrange(len(WatermarkedSubbandValue)):
        for j in xrange(len(WatermarkedSubbandValue[i])):
            WatermarkedSubbandValue[i][j] += OriginalSubbandValue[i][j]

    return WatermarkedSubbandValue

def power(A):
    summ = 0.0
    for i in A:
        summ += sum(i) ** 2
    return summ / (len(A) * len(A[0]))

def PSNR(img, img_):
    return 10 * math.log10(250**2 * (len(img) * len(img[0])) / sum(sum((img - img_)**2)))

def MSE(img, img_):
    return sum(sum((img - img_)**2)) / (len(img) * len(img[0]))

def calc_correlation_thresholds(LH, HL, HH, Pseudo):
    s = 0
    for i in xrange(len(Pseudo)):
        for j in xrange(len(Pseudo[i])):
            LH[i][j] *= Pseudo[i][j]
            HL[i][j] *= Pseudo[i][j]
            HH[i][j] *= Pseudo[i][j]

    for i in LH: s += sum(i)
    for i in HL: s += sum(i)
    for i in HH: s += sum(i)

    return s / (len(LH) * len(LH[0]) + len(HL) * len(HL[0]) + len(HH) * len(HH[0]))
    
# Prepare
image = Image.open('lena1.bmp').convert('L')
width, height = image.size
print width, height

wavelet = pywt.Wavelet('db1')      # use wavelet Dobeshi 4
image = numpy.array(image)
level = 4    # lwl decomposit
mode = 'sym' # protect mode
image_ = Image.open('/Lena/WMI_Gausse.png').convert('L')
image_ = numpy.array(image_)

# step 1
coeffs = pywt.wavedec2(image, wavelet, mode, level) # get coefs DWT
A, (LH4, HL4, HH4), (LH3, HL3, HH3), (LH2, HL2, HH2), (LH1, HL1, HH1) = coeffs

#print len(A), len(A[0]), len(LH4), len(LH4[0]), len(LH3), len(LH3[0]), len(LH2), len(LH2[0]), len(LH1), len(LH1[0])

# step 2
image_wm = Image.open('Watermark.bmp')#.convert('L') # load watermark
width_wm, height_wm = image_wm.size

# step 3
random_matrix = marking_rand_matr(image_wm, width_wm, height_wm)

# step 4
avgLH3, avgHL3, avgHH3, avgLH2, avgHL2, avgHH2 = 0, 0, 0, 0, 0, 0
avgA = 0

# find average value coeff subband
for i in LH3: avgLH3 += sum(i) / (len(i) * len(LH3))
for i in HL3: avgHL3 += sum(i) / (len(i) * len(HL3))
for i in HH3: avgHH3 += sum(i) / (len(i) * len(HH3))
for i in LH2: avgLH2 += sum(i) / (len(i) * len(LH2))
for i in HL2: avgHL2 += sum(i) / (len(i) * len(HL2))
for i in HH2: avgHH2 += sum(i) / (len(i) * len(HH2))
for i in A: avgA += sum(i) / (len(i) * len(A))

TotalChangeLH3 = [[0 for j in xrange(len(LH3[0]))] for i in xrange(len(LH3))]
TotalChangeHL3 = [[0 for j in xrange(len(HL3[0]))] for i in xrange(len(HL3))]
TotalChangeHH3 = [[0 for j in xrange(len(HH3[0]))] for i in xrange(len(HH3))]
TotalChangeLH2 = [[0 for j in xrange(len(LH2[0]))] for i in xrange(len(LH2))]
TotalChangeHL2 = [[0 for j in xrange(len(HL2[0]))] for i in xrange(len(HL2))]
TotalChangeHH2 = [[0 for j in xrange(len(HH2[0]))] for i in xrange(len(HH2))]

for i in xrange(len(A)):
    for j in xrange(len(A[i])):
        Brightness = avgA/1000
        if Brightness < 0.2: Brightness = 1 - Brightness
        
        if i == len(A) - 1 and j == len(A[i]) - 1:
            Var = varience(A[i][j], A[i][j], A[i][j], A[i][j], 4.0)
        elif i == len(A) - 1 and j != len(A[i]) - 1:
            Var = varience(A[i][j], A[i][j + 1], A[i][j], A[i][j + 1], 4.0)
        elif i != len(A) - 1 and j == len(A[i]) - 1:
            Var = varience(A[i][j], A[i + 1][j], A[i][j], A[i + 1][j], 4.0)
        else: Var = varience(A[i][j], A[i + 1][j],A[i][j + 1], A[i + 1][j + 1], 4.0)
        
        # # # # # #
        fg, L, a, b = True, 0.0625, 1, 1
        S = []
        for k in xrange(i * 2, i * 2 + 2):
            SS = []
            for l in xrange(j * 2, j * 2 + 2):
                if LH3[k][l] <= avgLH3: fg = False
                SS.append(LH3[k][l])
            S.append(SS)
            del(SS)

        if fg == True: Texture = L * Var * (power(S) / power(LH3))
        else: Texture = 1
        
        Texture = 1
        BlockChange = a * b * 0.16 * Brightness * Texture
        TotalChangeLH3[i][j] = BlockChange
        # # # # # #
        fg, L, a, b = True, 0.0625, 1, 1
        S = []
        for k in xrange(i * 2, i * 2 + 2):
            SS = []
            for l in xrange(j * 2, j * 2 + 2):
                if HL3[k][l] <= avgHL3: fg = False
                SS.append(HL3[k][l])
            S.append(SS)
            del(SS)
        
        if fg == True: Texture = L * Var * (power(S) / power(HL3))
        else: Texture = 1
        Texture = 1
        BlockChange = a * b * 0.16 * Brightness * Texture
        TotalChangeHL3[i][j] = BlockChange
        ######
        fg, L, a, b = True, 0.0625, 1, 2 ** 0.5
        S = []
        for k in xrange(i * 2, i * 2 + 2):
            SS = []
            for l in xrange(j * 2, j * 2 + 2):
                if HH3[k][l] <= avgHH3: fg = False
                SS.append(HH3[k][l])
            S.append(SS)
            del(SS)
        
        if fg == True: Texture = L * Var * (power(S) / power(HH3))
        else: Texture = 1
        Texture = 1
        BlockChange = a * b * 0.16 * Brightness * Texture
        TotalChangeHH3[i][j] = BlockChange
        # # # # # #
        fg, L, a, b = True, 1/256, 2, 1
        S = []
        for k in xrange(i * 4, i * 4 + 4):
            SS = []
            for l in xrange(j * 4, j * 4 + 4):
                if LH2[k][l] <= avgLH2: fg = False
                SS.append(LH2[k][l])
            S.append(SS)
            del(SS)
        
        if fg == True: Texture = L * Var * (power(S) / power(LH2))
        else: Texture = 1
        
        BlockChange = a * b * 0.16 * Brightness * Texture
        TotalChangeLH2[i][j] = BlockChange
        # # # # # #
        fg, L, a, b = True, 1/256, 2, 1
        S = []
        for k in xrange(i * 4, i * 4 + 4):
            SS = []
            for l in xrange(j * 4, j * 4 + 4):
                if HL2[k][l] <= avgHL2: fg = False
                SS.append(HL2[k][l])
            S.append(SS)
            del(SS)
        
        if fg == True: Texture = L * Var * (power(S) / power(HL2))
        else: Texture = 1
        
        BlockChange = a * b * 0.16 * Brightness * Texture
        TotalChangeHL2[i][j] = BlockChange
        # # # # # #
        fg, L, a, b = True, 1/256, 2, 2 ** 0.5
        S = []
        for k in xrange(i * 4, i * 4 + 4):
            SS = []
            for l in xrange(j * 4, j * 4 + 4):
                if HH2[k][l] <= avgHH2: fg = False
                SS.append(HH2[k][l])
            S.append(SS)
            del(SS)
        
        if fg == True: Texture = L * Var * (power(S) / power(HH2))
        else: Texture = 1
        
        BlockChange = a * b * 0.16 * Brightness * Texture
        TotalChangeHH2[i][j] = BlockChange
# step 5
w3 = 0.4
w2 = 0.4

WatermarkedLH3 = watermarked(LH3, w3, TotalChangeLH3, random_matrix)
WatermarkedHL3 = watermarked(HL3, w3, TotalChangeHL3, random_matrix)
WatermarkedHH3 = watermarked(HH3, w3, TotalChangeHH3, random_matrix)
WatermarkedLH2 = watermarked(LH2, w2, TotalChangeLH2, random_matrix)
WatermarkedHL2 = watermarked(HL2, w2, TotalChangeHL2, random_matrix)
WatermarkedHH2 = watermarked(HH2, w2, TotalChangeHH2, random_matrix)

lvlDWT = embeding_subbands(width, height, coeffs, level)

# step 6
Watermarked_coeffs = A, (LH4, HL4, HH4), (WatermarkedLH3, WatermarkedHL3, WatermarkedHH3),(WatermarkedLH2, WatermarkedHL2, WatermarkedHH2), (LH1, HL1, HH1)
ximage = pywt.waverec2(Watermarked_coeffs, wavelet, mode)

psnr = PSNR(image, ximage)
mse = MSE(image, ximage)
print 'psnr =', psnr
print 'mse =', mse

# step 7
T1 = calc_correlation_thresholds(LH3, HL3, HH3, random_matrix)
T2 = calc_correlation_thresholds(LH2, HL2, HH2, random_matrix)

coeffs = pywt.wavedec2(ximage, wavelet, mode, level) # get coefs DWT
#coeffs = pywt.wavedec2(image_, wavelet, mode, level)
A, (LH4, HL4, HH4), (LH3, HL3, HH3), (LH2, HL2, HH2), (LH1, HL1, HH1) = coeffs


# Detection WM
R1 = calc_correlation_thresholds(LH3, HL3, HH3, random_matrix)
R2 = calc_correlation_thresholds(LH2, HL2, HH2, random_matrix)

if R1 > T1 and R2 > T2: print 'Watermark detection'
elif R1 < T1 and R2 < T2: print 'Watermark detection'
else: print 'Watermark absent'


ii = Image.new('L', (len(ximage), len(ximage[0])))
for i in xrange(len(ximage)):
    for j in xrange(len(ximage[i])):
        ii.putpixel((j, i),ximage[i][j])

ii.save('WMI.png')

kk = ximage - image
plt.figure(1)
plt.imshow(ximage, cmap = plt.get_cmap('gray'))
plt.figure(2)
plt.imshow(kk, cmap = plt.get_cmap('gray'))
pylab.show()
