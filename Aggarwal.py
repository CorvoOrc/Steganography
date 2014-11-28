# -*- coding: cp1251 -*-

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
    for k in xrange(level):
        for i in xrange(W):
            data[H/2][i] = -1
        H, W = H/2, W/2
        
    W, H = w, h
    for k in xrange(level):
        for i in xrange(H):
            data[i][W/2] = -1
        H, W = H/2, W/2

    return data

def embeding_subbands(w, h, coeffs, level):
    LL1, (LH1, HL1, HH1) = coeffs
    matrixCoeffs = [[0 for j in xrange(w)] for i in xrange(h)]
    
    matrixCoeffs = embeding_subband(matrixCoeffs, 0,    w/2,  0,    h/2,  LL1)
    matrixCoeffs = embeding_subband(matrixCoeffs, 0,    w/2,  h/2,  h,    LH1)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/2,  w,    0,    h/2,  HL1)
    matrixCoeffs = embeding_subband(matrixCoeffs, w/2,  w,    h/2,  h,    HH1)

    matrixCoeffs = draw_edge(matrixCoeffs, w, h, level)
    
    return matrixCoeffs

def PSNR(img, img_):
    return 10 * math.log10(250**2 * (len(img) * len(img[0])) / sum(sum((img - img_)**2)))

def MSE(img, img_):
    return sum(sum((img - img_)**2)) / (len(img) * len(img[0]))

def SR(img, img_):
    S, D = 0.0, 0.0
    for i in xrange(len(img)):
        for j in xrange(len(img[i])):
            if int(img[i][j]) == int(img_[i][j]): S += 1
            else: D += 1

    return S / (S + D)

wavelet = pywt.Wavelet('db4')
level = 1    # lwl decomposition
mode = 'sym' # protected mode

name = raw_input('Input path to image:')
image = Image.open(name)
print 'Mode open image -', image.mode
if not image.mode in 'LIF': image = image.convert('L')
width, height = image.size
image = numpy.array(image)
#image_ = Image.open('Атаки/Peppers/WMI_noise.png').convert('L')
#image_ = numpy.array(image_)

name = raw_input('Enter primary WM:')
primaryWM = Image.open(name)
wp, hp = primaryWM.size
primaryWM = numpy.array(primaryWM)

coeffs = pywt.wavedec2(primaryWM, wavelet, mode, level)
cap1, (chp1, cvp1, cdp1) = coeffs
tree1 = embeding_subbands(wp, hp, coeffs, level)

name = raw_input('Enter secondary WM')
secondaryWM = Image.open(name)
secondaryWM = secondaryWM.resize((len(chp1), len(chp1[0])))
secondaryWM = numpy.array(secondaryWM)
fg = raw_input('Embeding at size or full(S-size/F-full)?')

cap1 += secondaryWM

coeffs = cap1, (chp1, cvp1, cdp1)
nestedWM = pywt.waverec2(coeffs, wavelet)

psnr = PSNR(primaryWM, nestedWM)
mse = MSE(primaryWM, nestedWM)
print 'nestedWM ws primaryWM'
print 'PSNR:', psnr, '\nMSE:', mse

coeffs = pywt.wavedec2(image, wavelet, mode, level)
ca1, (ch1, cv1, cd1) = coeffs
tree2 = embeding_subbands(width, height, coeffs, level)

if fg == 'S': nestedArray = numpy.array(nestedWM)
else:
    nestedww = Image.new('L', (len(nestedWM), len(nestedWM[0])))
    for i in xrange(len(nestedWM)):
        for j in xrange(len(nestedWM[i])):
            nestedww.putpixel((j, i), nestedWM[i][j])

    nestedww = nestedww.resize((len(ca1), len(ca1[0])))
    nestedArray = numpy.array(nestedww)

Ca1 = [[0 for j in xrange(len(ca1[0]))] for i in xrange(len(ca1))]
Cd1 = [[0 for j in xrange(len(cd1[0]))] for i in xrange(len(cd1))]

alp = 0.04
Ca1 += ca1
for i in xrange(len(nestedArray)):
    for j in xrange(len(nestedArray[i])):
        Ca1[i][j] += alp * nestedArray[i][j]

alp = 0.01
Cd1 += cd1
for i in xrange(len(nestedArray)):
    for j in xrange(len(nestedArray[i])):
        Cd1[i][j] += alp * nestedArray[i][j]

coeffs = Ca1, (ch1, cv1, Cd1)
WatermarkedImg = pywt.waverec2(coeffs, wavelet)

psnr = PSNR(image, WatermarkedImg)
mse = MSE(image, WatermarkedImg)
sr = SR(image, WatermarkedImg)
print '\nWatermarkedImg ws OriginalImg'
print 'PSNR:', psnr, '\nMSE:', mse, '\nSR:', sr

ii = Image.new('L', (len(WatermarkedImg), len(WatermarkedImg[0])))
for i in xrange(len(WatermarkedImg)):
    for j in xrange(len(WatermarkedImg[i])):
        ii.putpixel((j, i),WatermarkedImg[i][j])

ii.save('WMI.png')

###### Watermark Extraction
coeffs = pywt.wavedec2(WatermarkedImg, wavelet, mode, level)
#coeffs = pywt.wavedec2(image_, wavelet, mode, level)
c1a1, (c1h1, c1v1, c1d1) = coeffs

LLw = [[0 for j in xrange(len(ca1[0]))] for i in xrange(len(ca1))]
HHw = [[0 for j in xrange(len(cd1[0]))] for i in xrange(len(cd1))]

alp = 0.04
LLw = (c1a1 - ca1) / alp
    
alp = 0.01
HHw = (c1d1 - cd1) / alp

print "Emdeging success"

plt.figure(1)
plt.title('Level-1 DWT primary WM')
plt.imshow(tree1, cmap = plt.get_cmap('gray'))
plt.figure(2)
plt.title('Nested WM')
plt.imshow(nestedWM, cmap = plt.get_cmap('gray'))
plt.figure(3)
plt.title('Level-1 DWT cover image')
plt.imshow(tree2, cmap = plt.get_cmap('gray'))
plt.figure(4)
plt.title('Watermarked Image')
plt.imshow(WatermarkedImg, cmap = plt.get_cmap('gray'))
plt.figure(5)
plt.title('Extracting WM in LL band')
plt.imshow(LLw, cmap = plt.get_cmap('gray'))
plt.figure(6)
plt.title('Extracting WM in HH band')
plt.imshow(HHw, cmap = plt.get_cmap('gray'))
pylab.show()
