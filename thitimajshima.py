__author__ = 'Steshenko'

# -*- coding: cp1251 -*-

import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
import matplotlib.pyplot as plt
import matplotlib
import pywt, numpy
import pylab
import math
import random


def embeding_subband(data, w_b, w_e, h_b, h_e, coef):
    k, l = 0, 0
    for i in xrange(h_b, h_e):
        for j in xrange(w_b, w_e):
            data[i][j], l = coef[k][l], l + 1
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
    s, d = 0.0, 0.0
    for i in xrange(len(img)):
        for j in xrange(len(img[i])):
            if int(img[i][j]) == int(img_[i][j]):
                s += 1
            else: d += 1

    return s / (s + d)


def pdf(x):
    return (1 / (math.sqrt(2*math.pi))) * math.exp(-(x**2) / 2)

wavelet = pywt.Wavelet('db4')
level = 2
mode = 'sym'

name = 'hill.png'  # raw_input('Enter path to image:')
image = Image.open(name)
print 'Mode open image -', image.mode
if not (image.mode in 'LIF'):
    image = image.convert('L')
width, height = image.size
image = numpy.array(image)
image_ = Image.open('cat.jpg').convert('L')
image_ = numpy.array(image_)

coeffs = pywt.wavedec2(image, wavelet, mode, level)
ca2, (ch2, cv2, cd2), (ch1, cv1, cd1) = coeffs

C1, C2 = 0.0, 0.0

for i in ch1:
    for j in i:
        if j > C1:
            C1 = j
for i in cv1:
    for j in i:
        if j > C1:
            C1 = j
for i in cd1:
    for j in i:
        if j > C1:
            C1 = j

for i in ch2:
    for j in i:
        if j > C2:
            C2 = j
for i in cv2:
    for j in i:
        if j > C2:
            C2 = j
for i in cd2:
    for j in i:
        if j > C2:
            C2 = j

T1 = 2 ** (int(math.log(C1, 2)) - 1)
T2 = 2 ** (int(math.log(C2, 2)) - 1)
print 'T1={t1} T2={t2} C1={c1} C2={c2}'.format(t1=T1, t2=T2, c1=C1, c2=C2)
PRS = []
for i in xrange(2000):
    PRS.append(pdf(random.gauss(0, 1)))

alp, x = 0.2, 0
for i in xrange(len(ch1)):
    for j in xrange(len(ch1[i])):
        if ch1[i][j] > T1:
            ch1[i][j] += alp * abs(ch1[i][j]) * PRS[x]
            x += 1

alp, x = 0.2, 0
for i in xrange(len(cv1)):
    for j in xrange(len(cv1[i])):
        if cv1[i][j] > T1:
            cv1[i][j] += alp * abs(cv1[i][j]) * PRS[x]
            x += 1

alp, x = 0.2, 0
for i in xrange(len(cd1)):
    for j in xrange(len(cd1[i])):
        if cd1[i][j] > T1:
            cd1[i][j] += alp * abs(cd1[i][j]) * PRS[x]
            x += 1

alp, x = 0.2, 0
for i in xrange(len(ch2)):
    for j in xrange(len(ch2[i])):
        if ch2[i][j] > T2:
            ch2[i][j] += alp * abs(ch2[i][j]) * PRS[x]
            x += 1

alp, x = 0.2, 0
for i in xrange(len(cv2)):
    for j in xrange(len(cv2[i])):
        if cv2[i][j] > T2:
            cv2[i][j] += alp * abs(cv2[i][j]) * PRS[x]
            x += 1

alp, x = 0.2, 0
for i in xrange(len(cd2)):
    for j in xrange(len(cd2[i])):
        if cd2[i][j] > T2:
            cd2[i][j] += alp * abs(cd2[i][j]) * PRS[x]
            x += 1

coeffs = ca2, (ch2, cv2, cd2), (ch1, cv1, cd1)
WatermarkedImg = pywt.waverec2(coeffs, wavelet)

iris = WatermarkedImg - image
psnr = PSNR(image, WatermarkedImg)
mse = MSE(image, WatermarkedImg)
print 'PSNR={psnr}\nMSE={mse}'.format(psnr=psnr, mse=mse)

ii = Image.new('L', (len(WatermarkedImg), len(WatermarkedImg[0])))
for i in xrange(len(WatermarkedImg)):
    for j in xrange(len(WatermarkedImg[i])):
        ii.putpixel((j, i),WatermarkedImg[i][j])

ii.save('WMI.png')

# Watermark Detection
coeffs = pywt.wavedec2(WatermarkedImg, wavelet, mode, level)
ca2, (ch2, cv2, cd2), (ch1, cv1, cd1) = coeffs

C1, C2 = 0.0, 0.0
for i in ch1:
    for j in i:
        if j > C1:
            C1 = j
for i in cv1:
    for j in i:
        if j > C1:
            C1 = j
for i in cd1:
    for j in i:
        if j > C1:
            C1 = j

for i in ch2:
    for j in i:
        if j > C2:
            C2 = j
for i in cv2:
    for j in i:
        if j > C2:
            C2 = j
for i in cd2:
    for j in i:
        if j > C2:
            C2 = j

T1 = 2 ** (int(math.log(C1, 2)) - 1)
T2 = 2 ** (int(math.log(C2, 2)) - 1)

print 'T1={t1} T2={t2} C1={c1} C2={c2}'.format(t1=T1, t2=T2, c1=C1, c2=C2)
PRS = []
for i in xrange(2000):
    PRS.append(pdf(random.gauss(0, 1)))

z, S = 0.0, 0.0
alp, y = 0.2, 0
M = 0
for i in xrange(len(ch1)):
    for j in xrange(len(ch1[i])):
        if ch1[i][j] > T1:
            z += abs(ch1[i][j]) * PRS[y]
            S += abs(ch1[i][j])
            y += 1
M += y
y = 0
for i in xrange(len(cv1)):
    for j in xrange(len(cv1[i])):
        if cv1[i][j] > T1:
            z += abs(cv1[i][j]) * PRS[y]
            S += abs(cv1[i][j])
            y += 1

M += y
y = 0
for i in xrange(len(cd1)):
    for j in xrange(len(cd1[i])):
        if cd1[i][j] > T1:
            z += abs(cd1[i][j]) * PRS[y]
            S += abs(cd1[i][j])
            y += 1

M += y
y = 0
for i in xrange(len(ch2)):
    for j in xrange(len(ch2[i])):
        if ch2[i][j] > T2:
            z += abs(ch2[i][j]) * PRS[y]
            S += abs(ch2[i][j])
            y += 1

M += y
y = 0
for i in xrange(len(cv2)):
    for j in xrange(len(cv2[i])):
        if cv2[i][j] > T2:
            z += abs(cv2[i][j]) * PRS[y]
            S += abs(cv2[i][j])
            y += 1

M += y
y = 0
for i in xrange(len(cd2)):
    for j in xrange(len(cd2[i])):
        if cd2[i][j] > T2:
            z += abs(cd2[i][j]) * PRS[y]
            S += abs(cd2[i][j])
            y += 1
M += y
z /= M
S *= (alp / (2*M))

print 'z={z} S={S}'.format(z=z, S=S)
if z > S:
    print 'Watermark present'
else:
    print 'Watermark absent'

plt.figure(1)
plt.title('Original image')
plt.imshow(image, cmap=plt.get_cmap('gray'))
plt.figure(2)
plt.title('Watermarked image')
plt.imshow(WatermarkedImg, cmap=plt.get_cmap('gray'))
plt.figure(3)
plt.title('Subtract')
plt.imshow(iris, cmap=plt.get_cmap('gray'))
pylab.show()
