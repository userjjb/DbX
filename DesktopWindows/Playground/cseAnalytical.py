import numpy as np
#import sympy
#import sympy.abc
#
#atan = sympy.atan
#log = sympy.log
#sqrt = sympy.sqrt
#pi = sympy.pi
#I =sympy.I
#
#a,b,c,d,e,f = sympy.symbols('a b c d e f ')
#repl, redu = sympy.cse([(4*e**2*atan((a*c)/(e*sqrt(a**2 + c**2 + e**2))) - 4*e**2*atan((b*c)/(e*sqrt(b**2 + c**2 + e**2))) - 4*e**2*atan((a*d)/(e*sqrt(a**2 + d**2 + e**2))) + 4*e**2*atan((b*d)/(e*sqrt(b**2 + d**2 + e**2))) + 2*(2*c**2*(atan((a*e)/(c*sqrt(a**2 + c**2 + e**2))) - atan((b*e)/(c*sqrt(b**2 + c**2 + e**2))) - atan((a*f)/(c*sqrt(a**2 + c**2 + f**2))) + atan((b*f)/(c*sqrt(b**2 + c**2 + f**2)))) + a**2*(atan((c*e)/(a*sqrt(a**2 + c**2 + e**2))) - atan((d*e)/(a*sqrt(a**2 + d**2 + e**2))) - atan((c*f)/(a*sqrt(a**2 + c**2 + f**2))) + atan((d*f)/(a*sqrt(a**2 + d**2 + f**2)))) + 2*(2*(c**2 + e**2)*pi + f**2*(-atan((a*c)/(f*sqrt(a**2 + c**2 + f**2))) + atan((b*c)/(f*sqrt(b**2 + c**2 + f**2))) + atan((a*d)/(f*sqrt(a**2 + d**2 + f**2))) - atan((b*d)/(f*sqrt(b**2 + d**2 + f**2)))) + d**2*(-atan((a*e)/(d*sqrt(a**2 + d**2 + e**2))) + atan((b*e)/(d*sqrt(b**2 + d**2 + e**2))) + atan((a*f)/(d*sqrt(a**2 + d**2 + f**2))) - atan((b*f)/(d*sqrt(b**2 + d**2 + f**2))))) + b**2*(-atan((c*e)/(b*sqrt(b**2 + c**2 + e**2))) + atan((d*e)/(b*sqrt(b**2 + d**2 + e**2))) + atan((c*f)/(b*sqrt(b**2 + c**2 + f**2))) - atan((d*f)/(b*sqrt(b**2 + d**2 + f**2))))) - 4*a*e*log(a**2 + e**2) + 4*b*e*log(b**2 + e**2) - 4*c*e*log(c**2 + e**2) + 4*d*e*log(d**2 + e**2) + 4*c*e*log(-a + sqrt(a**2 + c**2 + e**2)) + 4*a*e*log(-c + sqrt(a**2 + c**2 + e**2)) + 4*c*e*log(b + sqrt(b**2 + c**2 + e**2)) - 4*b*e*log(-c + sqrt(b**2 + c**2 + e**2)) - 4*d*e*log(-a + sqrt(a**2 + d**2 + e**2)) + 4*a*e*log(d + sqrt(a**2 + d**2 + e**2)) - 4*d*e*log(b + sqrt(b**2 + d**2 + e**2)) - 4*b*e*log(d + sqrt(b**2 + d**2 + e**2)) - I*e**2*log(((a - I*e)*e - I*c*(c + sqrt(a**2 + c**2 + e**2)))/(a - I*e)) + I*e**2*log(((a + I*e)*e + I*c*(c + sqrt(a**2 + c**2 + e**2)))/(a + I*e)) - I*c**2*log((I*a*c + c**2 + e*(e + sqrt(a**2 + c**2 + e**2)))/(I*a + c)) + I*e**2*log(((b - I*e)*e - I*c*(c + sqrt(b**2 + c**2 + e**2)))/(b - I*e)) - I*e**2*log(((b + I*e)*e + I*c*(c + sqrt(b**2 + c**2 + e**2)))/(b + I*e)) + I*c**2*log((I*b*c + c**2 + e*(e + sqrt(b**2 + c**2 + e**2)))/(I*b + c)) + I*e**2*log(((a - I*e)*e - I*d*(d + sqrt(a**2 + d**2 + e**2)))/(a - I*e)) - I*e**2*log(((a + I*e)*e + I*d*(d + sqrt(a**2 + d**2 + e**2)))/(a + I*e)) + I*d**2*log((I*a*d + d**2 + e*(e + sqrt(a**2 + d**2 + e**2)))/(I*a + d)) - I*e**2*log(((b - I*e)*e - I*d*(d + sqrt(b**2 + d**2 + e**2)))/(b - I*e)) + I*e**2*log(((b + I*e)*e + I*d*(d + sqrt(b**2 + d**2 + e**2)))/(b + I*e)) - I*d**2*log((I*b*d + d**2 + e*(e + sqrt(b**2 + d**2 + e**2)))/(I*b + d)) + I*c**2*log((a*c + I*(c**2 + e*(e + sqrt(a**2 + c**2 + e**2))))/(a + I*c)) - I*c**2*log((b*c + I*(c**2 + e*(e + sqrt(b**2 + c**2 + e**2))))/(b + I*c)) - I*d**2*log((a*d + I*(d**2 + e*(e + sqrt(a**2 + d**2 + e**2))))/(a + I*d)) + I*d**2*log((b*d + I*(d**2 + e*(e + sqrt(b**2 + d**2 + e**2))))/(b + I*d)) + 4*a*f*log(a**2 + f**2) - 4*b*f*log(b**2 + f**2) + 4*c*f*log(c**2 + f**2) - 4*d*f*log(d**2 + f**2) - 4*c*f*log(-a + sqrt(a**2 + c**2 + f**2)) - 4*a*f*log(-c + sqrt(a**2 + c**2 + f**2)) + 4*a*c*log((f + sqrt(a**2 + c**2 + f**2))/(e + sqrt(a**2 + c**2 + e**2))) - 4*c*f*log(b + sqrt(b**2 + c**2 + f**2)) + 4*b*f*log(-c + sqrt(b**2 + c**2 + f**2)) - 4*b*c*log((f + sqrt(b**2 + c**2 + f**2))/(e + sqrt(b**2 + c**2 + e**2))) + 4*d*f*log(-a + sqrt(a**2 + d**2 + f**2)) - 4*a*f*log(d + sqrt(a**2 + d**2 + f**2)) - 4*a*d*log((f + sqrt(a**2 + d**2 + f**2))/(e + sqrt(a**2 + d**2 + e**2))) + 4*d*f*log(b + sqrt(b**2 + d**2 + f**2)) + 4*b*f*log(d + sqrt(b**2 + d**2 + f**2)) + 4*b*d*log((f + sqrt(b**2 + d**2 + f**2))/(e + sqrt(b**2 + d**2 + e**2))) + I*f**2*log(((a - I*f)*f - I*c*(c + sqrt(a**2 + c**2 + f**2)))/(a - I*f)) - I*f**2*log(((a + I*f)*f + I*c*(c + sqrt(a**2 + c**2 + f**2)))/(a + I*f)) + I*c**2*log((I*a*c + c**2 + f*(f + sqrt(a**2 + c**2 + f**2)))/(I*a + c)) - I*f**2*log(((b - I*f)*f - I*c*(c + sqrt(b**2 + c**2 + f**2)))/(b - I*f)) + I*f**2*log(((b + I*f)*f + I*c*(c + sqrt(b**2 + c**2 + f**2)))/(b + I*f)) - I*c**2*log((I*b*c + c**2 + f*(f + sqrt(b**2 + c**2 + f**2)))/(I*b + c)) - I*f**2*log(((a - I*f)*f - I*d*(d + sqrt(a**2 + d**2 + f**2)))/(a - I*f)) + I*f**2*log(((a + I*f)*f + I*d*(d + sqrt(a**2 + d**2 + f**2)))/(a + I*f)) - I*d**2*log((I*a*d + d**2 + f*(f + sqrt(a**2 + d**2 + f**2)))/(I*a + d)) + I*f**2*log(((b - I*f)*f - I*d*(d + sqrt(b**2 + d**2 + f**2)))/(b - I*f)) - I*f**2*log(((b + I*f)*f + I*d*(d + sqrt(b**2 + d**2 + f**2)))/(b + I*f)) + I*d**2*log((I*b*d + d**2 + f*(f + sqrt(b**2 + d**2 + f**2)))/(I*b + d)) - I*c**2*log((a*c + I*(c**2 + f*(f + sqrt(a**2 + c**2 + f**2))))/(a + I*c)) + I*c**2*log((b*c + I*(c**2 + f*(f + sqrt(b**2 + c**2 + f**2))))/(b + I*c)) + I*d**2*log((a*d + I*(d**2 + f*(f + sqrt(a**2 + d**2 + f**2))))/(a + I*d)) - I*d**2*log((b*d + I*(d**2 + f*(f + sqrt(b**2 + d**2 + f**2))))/(b + I*d)))/4])
#
#repl2, redu2 = sympy.cse([(4*e**2*atan((a*c)/(e*sqrt(a**2 + c**2 + e**2))) - 4*e**2*atan((b*c)/(e*sqrt(b**2 + c**2 + e**2))) - 4*e**2*atan((a*d)/(e*sqrt(a**2 + d**2 + e**2))) + 4*e**2*atan((b*d)/(e*sqrt(b**2 + d**2 + e**2))) + 2*(2*c**2*(atan((a*e)/(c*sqrt(a**2 + c**2 + e**2))) - atan((b*e)/(c*sqrt(b**2 + c**2 + e**2))) - atan((a*f)/(c*sqrt(a**2 + c**2 + f**2))) + atan((b*f)/(c*sqrt(b**2 + c**2 + f**2)))) + a**2*(atan((c*e)/(a*sqrt(a**2 + c**2 + e**2))) - atan((d*e)/(a*sqrt(a**2 + d**2 + e**2))) - atan((c*f)/(a*sqrt(a**2 + c**2 + f**2))) + atan((d*f)/(a*sqrt(a**2 + d**2 + f**2)))) + 2*(2*(c**2 + e**2)*pi + f**2*(-atan((a*c)/(f*sqrt(a**2 + c**2 + f**2))) + atan((b*c)/(f*sqrt(b**2 + c**2 + f**2))) + atan((a*d)/(f*sqrt(a**2 + d**2 + f**2))) - atan((b*d)/(f*sqrt(b**2 + d**2 + f**2)))) + d**2*(-atan((a*e)/(d*sqrt(a**2 + d**2 + e**2))) + atan((b*e)/(d*sqrt(b**2 + d**2 + e**2))) + atan((a*f)/(d*sqrt(a**2 + d**2 + f**2))) - atan((b*f)/(d*sqrt(b**2 + d**2 + f**2))))) + b**2*(-atan((c*e)/(b*sqrt(b**2 + c**2 + e**2))) + atan((d*e)/(b*sqrt(b**2 + d**2 + e**2))) + atan((c*f)/(b*sqrt(b**2 + c**2 + f**2))) - atan((d*f)/(b*sqrt(b**2 + d**2 + f**2))))) - 4*a*e*log(a**2 + e**2) + 4*b*e*log(b**2 + e**2) - 4*c*e*log(c**2 + e**2) + 4*d*e*log(d**2 + e**2) + 4*c*e*log(-a + sqrt(a**2 + c**2 + e**2)) + 4*a*e*log(-c + sqrt(a**2 + c**2 + e**2)) + 4*c*e*log(b + sqrt(b**2 + c**2 + e**2)) - 4*b*e*log(-c + sqrt(b**2 + c**2 + e**2)) - 4*d*e*log(-a + sqrt(a**2 + d**2 + e**2)) + 4*a*e*log(d + sqrt(a**2 + d**2 + e**2)) - 4*d*e*log(b + sqrt(b**2 + d**2 + e**2)) - 4*b*e*log(d + sqrt(b**2 + d**2 + e**2)) - I*e**2*log(((a - I*e)*e - I*c*(c + sqrt(a**2 + c**2 + e**2)))/(a - I*e)) + I*e**2*log(((a + I*e)*e + I*c*(c + sqrt(a**2 + c**2 + e**2)))/(a + I*e)) - I*c**2*log((I*a*c + c**2 + e*(e + sqrt(a**2 + c**2 + e**2)))/(I*a + c)) + I*e**2*log(((b - I*e)*e - I*c*(c + sqrt(b**2 + c**2 + e**2)))/(b - I*e)) - I*e**2*log(((b + I*e)*e + I*c*(c + sqrt(b**2 + c**2 + e**2)))/(b + I*e)) + I*c**2*log((I*b*c + c**2 + e*(e + sqrt(b**2 + c**2 + e**2)))/(I*b + c)) + I*e**2*log(((a - I*e)*e - I*d*(d + sqrt(a**2 + d**2 + e**2)))/(a - I*e)) - I*e**2*log(((a + I*e)*e + I*d*(d + sqrt(a**2 + d**2 + e**2)))/(a + I*e)) + I*d**2*log((I*a*d + d**2 + e*(e + sqrt(a**2 + d**2 + e**2)))/(I*a + d)) - I*e**2*log(((b - I*e)*e - I*d*(d + sqrt(b**2 + d**2 + e**2)))/(b - I*e)) + I*e**2*log(((b + I*e)*e + I*d*(d + sqrt(b**2 + d**2 + e**2)))/(b + I*e)) - I*d**2*log((I*b*d + d**2 + e*(e + sqrt(b**2 + d**2 + e**2)))/(I*b + d)) + I*c**2*log((a*c + I*(c**2 + e*(e + sqrt(a**2 + c**2 + e**2))))/(a + I*c)) - I*c**2*log((b*c + I*(c**2 + e*(e + sqrt(b**2 + c**2 + e**2))))/(b + I*c)) - I*d**2*log((a*d + I*(d**2 + e*(e + sqrt(a**2 + d**2 + e**2))))/(a + I*d)) + I*d**2*log((b*d + I*(d**2 + e*(e + sqrt(b**2 + d**2 + e**2))))/(b + I*d)) + 4*a*f*log(a**2 + f**2) - 4*b*f*log(b**2 + f**2) + 4*c*f*log(c**2 + f**2) - 4*d*f*log(d**2 + f**2) - 4*c*f*log(-a + sqrt(a**2 + c**2 + f**2)) - 4*a*f*log(-c + sqrt(a**2 + c**2 + f**2)) + 4*a*c*log((f + sqrt(a**2 + c**2 + f**2))/(e + sqrt(a**2 + c**2 + e**2))) - 4*c*f*log(b + sqrt(b**2 + c**2 + f**2)) + 4*b*f*log(-c + sqrt(b**2 + c**2 + f**2)) - 4*b*c*log((f + sqrt(b**2 + c**2 + f**2))/(e + sqrt(b**2 + c**2 + e**2))) + 4*d*f*log(-a + sqrt(a**2 + d**2 + f**2)) - 4*a*f*log(d + sqrt(a**2 + d**2 + f**2)) - 4*a*d*log((f + sqrt(a**2 + d**2 + f**2))/(e + sqrt(a**2 + d**2 + e**2))) + 4*d*f*log(b + sqrt(b**2 + d**2 + f**2)) + 4*b*f*log(d + sqrt(b**2 + d**2 + f**2)) + 4*b*d*log((f + sqrt(b**2 + d**2 + f**2))/(e + sqrt(b**2 + d**2 + e**2))) + I*f**2*log(((a - I*f)*f - I*c*(c + sqrt(a**2 + c**2 + f**2)))/(a - I*f)) - I*f**2*log(((a + I*f)*f + I*c*(c + sqrt(a**2 + c**2 + f**2)))/(a + I*f)) + I*c**2*log((I*a*c + c**2 + f*(f + sqrt(a**2 + c**2 + f**2)))/(I*a + c)) - I*f**2*log(((b - I*f)*f - I*c*(c + sqrt(b**2 + c**2 + f**2)))/(b - I*f)) + I*f**2*log(((b + I*f)*f + I*c*(c + sqrt(b**2 + c**2 + f**2)))/(b + I*f)) - I*c**2*log((I*b*c + c**2 + f*(f + sqrt(b**2 + c**2 + f**2)))/(I*b + c)) - I*f**2*log(((a - I*f)*f - I*d*(d + sqrt(a**2 + d**2 + f**2)))/(a - I*f)) + I*f**2*log(((a + I*f)*f + I*d*(d + sqrt(a**2 + d**2 + f**2)))/(a + I*f)) - I*d**2*log((I*a*d + d**2 + f*(f + sqrt(a**2 + d**2 + f**2)))/(I*a + d)) + I*f**2*log(((b - I*f)*f - I*d*(d + sqrt(b**2 + d**2 + f**2)))/(b - I*f)) - I*f**2*log(((b + I*f)*f + I*d*(d + sqrt(b**2 + d**2 + f**2)))/(b + I*f)) + I*d**2*log((I*b*d + d**2 + f*(f + sqrt(b**2 + d**2 + f**2)))/(I*b + d)) - I*c**2*log((a*c + I*(c**2 + f*(f + sqrt(a**2 + c**2 + f**2))))/(a + I*c)) + I*c**2*log((b*c + I*(c**2 + f*(f + sqrt(b**2 + c**2 + f**2))))/(b + I*c)) + I*d**2*log((a*d + I*(d**2 + f*(f + sqrt(a**2 + d**2 + f**2))))/(a + I*d)) - I*d**2*log((b*d + I*(d**2 + f*(f + sqrt(b**2 + d**2 + f**2))))/(b + I*d)))/4], optimizations = 'basic')
#
#funs = []
#syms = [a,b,c,d,e,f]
#for i, v in enumerate(repl):
#    funs.append(sp.lambdify(syms,v[1],modules='numpy'))
#    syms.append(v[0])
#    
#glam = sp.lambdify(syms,redu[0],modules='numpy')

def exact_u(bounds):
    a, b, c, d, e, f = bounds
    sqrt = np.sqrt
    log = np.log
    atan = np.arctan
    pi = np.pi

    x0 = a*f
    x1 = a**2
    x2 = f**2
    x3 = x1 + x2
    x4 = b*e
    x5 = b**2
    x6 = e**2
    x7 = x5 + x6
    x8 = c*f
    x9 = c**2
    x10 = d*e
    x11 = d**2
    x12 = a*e
    x13 = x1 + x6
    x14 = b*f
    x15 = x2 + x5
    x16 = c*e
    x17 = d*f
    x18 = sqrt(x11 + x13)
    x19 = d + x18
    x20 = sqrt(x11 + x15)
    x21 = d + x20
    x22 = sqrt(x7 + x9)
    x23 = sqrt(x11 + x3)
    x24 = d + x23
    x25 = sqrt(x11 + x7)
    x26 = d + x25
    x27 = sqrt(x15 + x9)
    x28 = -c
    x29 = sqrt(x13 + x9)
    x30 = -a
    x31 = sqrt(x3 + x9)
    x32 = a*c
    x33 = 1/e
    x34 = 1/x29
    x35 = b*d
    x36 = 1/x25
    x37 = a*d
    x38 = 1/x18
    x39 = b*c
    x40 = 1/x22
    x41 = f + x31
    x42 = e + x29
    x43 = f + x20
    x44 = e + x25
    x45 = f + x23
    x46 = e + x18
    x47 = f + x27
    x48 = e + x22
    x49 = 1j*x9/4
    x50 = 1j*a
    x51 = 1/(c + x50)
    x52 = c*x50
    x53 = e*x42 + x9
    x54 = f*x41 + x9
    x55 = 1j*b
    x56 = 1/(c + x55)
    x57 = c*x55
    x58 = e*x48 + x9
    x59 = f*x47 + x9
    x60 = 1j*x11/4
    x61 = 1/(d + x50)
    x62 = d*x50
    x63 = e*x46 + x11
    x64 = f*x45 + x11
    x65 = 1/(d + x55)
    x66 = d*x55
    x67 = e*x44 + x11
    x68 = f*x43 + x11
    x69 = 1j*x6/4
    x70 = 1j*e
    x71 = a + x70
    x72 = 1/x71
    x73 = e*x71
    x74 = 1j*c
    x75 = x74*(c + x29)
    x76 = 1j*d
    x77 = x19*x76
    x78 = b + x70
    x79 = 1/x78
    x80 = e*x78
    x81 = x74*(c + x22)
    x82 = x26*x76
    x83 = 1j*x2/4
    x84 = 1j*f
    x85 = a + x84
    x86 = 1/x85
    x87 = f*x85
    x88 = x74*(c + x31)
    x89 = x24*x76
    x90 = b + x84
    x91 = 1/x90
    x92 = f*x90
    x93 = x74*(c + x27)
    x94 = x21*x76
    x95 = 1/(a + x74)
    x96 = 1/(b + x74)
    x97 = 1/(a + x76)
    x98 = 1/(b + x76)
    x99 = -x70
    x100 = a + x99
    x101 = 1/x100
    x102 = e*x100
    x103 = b + x99
    x104 = 1/x103
    x105 = e*x103
    x106 = -x84
    x107 = a + x106
    x108 = 1/x107
    x109 = f*x107
    x110 = b + x106
    x111 = 1/x110
    x112 = f*x110
    x113 = 1/c
    x114 = 1/x27
    x115 = 1/x31
    x116 = 1/d
    x117 = 1/x23
    x118 = 1/x20
    x119 = 1/f
    x120 = 1/a
    x121 = 1/b
    
    return np.real(-x0*log(x24) + x0*log(x3) - x0*log(x28 + x31) + x1*(-atan(x10*x120*x38) - atan(x115*x120*x8) + atan(x117*x120*x17) + atan(x120*x16*x34))/2 - x10*log(b + x25) + x10*log(x11 + x6) - x10*log(x18 + x30) + x11*(atan(x0*x116*x117) - atan(x116*x118*x14) - atan(x116*x12*x38) + atan(x116*x36*x4)) - x12*log(x13) + x12*log(x19) + x12*log(x28 + x29) - x14*log(x15) + x14*log(x21) + x14*log(x27 + x28) + x16*log(b + x22) + x16*log(x29 + x30) - x16*log(x6 + x9) + x17*log(b + x20) - x17*log(x11 + x2) + x17*log(x23 + x30) + x2*(atan(x114*x119*x39) - atan(x115*x119*x32) + atan(x117*x119*x37) - atan(x118*x119*x35)) + x32*log(x41/x42) + x35*log(x43/x44) - x37*log(x45/x46) - x39*log(x47/x48) - x4*log(x26) + x4*log(x7) - x4*log(x22 + x28) - x49*log(x51*(x52 + x53)) + x49*log(x51*(x52 + x54)) + x49*log(x56*(x57 + x58)) - x49*log(x56*(x57 + x59)) + x49*log(x95*(x32 + 1j*x53)) - x49*log(x95*(x32 + 1j*x54)) - x49*log(x96*(x39 + 1j*x58)) + x49*log(x96*(x39 + 1j*x59)) + x5*(atan(x10*x121*x36) + atan(x114*x121*x8) - atan(x118*x121*x17) - atan(x121*x16*x40))/2 + x6*atan(x32*x33*x34) + x6*atan(x33*x35*x36) - x6*atan(x33*x37*x38) - x6*atan(x33*x39*x40) + x60*log(x61*(x62 + x63)) - x60*log(x61*(x62 + x64)) - x60*log(x65*(x66 + x67)) + x60*log(x65*(x66 + x68)) - x60*log(x97*(x37 + 1j*x63)) + x60*log(x97*(x37 + 1j*x64)) + x60*log(x98*(x35 + 1j*x67)) - x60*log(x98*(x35 + 1j*x68)) - x69*log(x101*(x102 - x75)) + x69*log(x101*(x102 - x77)) + x69*log(x104*(x105 - x81)) - x69*log(x104*(x105 - x82)) + x69*log(x72*(x73 + x75)) - x69*log(x72*(x73 + x77)) - x69*log(x79*(x80 + x81)) + x69*log(x79*(x80 + x82)) - x8*log(b + x27) + x8*log(x2 + x9) - x8*log(x30 + x31) + x83*log(x108*(x109 - x88)) - x83*log(x108*(x109 - x89)) - x83*log(x111*(x112 - x93)) + x83*log(x111*(x112 - x94)) - x83*log(x86*(x87 + x88)) + x83*log(x86*(x87 + x89)) + x83*log(x91*(x92 + x93)) - x83*log(x91*(x92 + x94)) + x9*(-atan(x0*x113*x115) + atan(x113*x114*x14) + atan(x113*x12*x34) - atan(x113*x4*x40)) + pi*(2*x6 + 2*x9))

def exact_u2(bounds):
    a, b, c, d, e, f = bounds
    sqrt = np.sqrt
    log = np.log
    atan = np.arctan
    pi = np.pi

    x0 = c**2
    x1 = e**2
    x2 = x0 + x1
    x3 = a*f
    x4 = a**2
    x5 = f**2
    x6 = x4 + x5
    x7 = b*e
    x8 = b**2
    x9 = x1 + x8
    x10 = c*f
    x11 = d*e
    x12 = d**2
    x13 = a*e
    x14 = x1 + x4
    x15 = b*f
    x16 = x5 + x8
    x17 = c*e
    x18 = d*f
    x19 = sqrt(x12 + x14)
    x20 = d + x19
    x21 = sqrt(x12 + x16)
    x22 = d + x21
    x23 = sqrt(x0 + x9)
    x24 = sqrt(x12 + x6)
    x25 = d + x24
    x26 = sqrt(x12 + x9)
    x27 = d + x26
    x28 = sqrt(x0 + x16)
    x29 = -c
    x30 = sqrt(x0 + x14)
    x31 = -a
    x32 = sqrt(x0 + x6)
    x33 = a*c
    x34 = 1/e
    x35 = 1/x30
    x36 = b*d
    x37 = 1/x26
    x38 = a*d
    x39 = 1/x19
    x40 = b*c
    x41 = 1/x23
    x42 = f + x32
    x43 = e + x30
    x44 = f + x21
    x45 = e + x26
    x46 = f + x24
    x47 = e + x19
    x48 = f + x28
    x49 = e + x23
    x50 = 1j*x0/4
    x51 = 1j*a
    x52 = 1/(c + x51)
    x53 = c*x51
    x54 = e*x43 + x0
    x55 = f*x42 + x0
    x56 = 1j*b
    x57 = 1/(c + x56)
    x58 = c*x56
    x59 = e*x49 + x0
    x60 = f*x48 + x0
    x61 = 1j*x12/4
    x62 = 1/(d + x51)
    x63 = d*x51
    x64 = e*x47 + x12
    x65 = f*x46 + x12
    x66 = 1/(d + x56)
    x67 = d*x56
    x68 = e*x45 + x12
    x69 = f*x44 + x12
    x70 = 1j*x1/4
    x71 = 1j*e
    x72 = a + x71
    x73 = 1/x72
    x74 = e*x72
    x75 = 1j*c
    x76 = x75*(c + x30)
    x77 = 1j*d
    x78 = x20*x77
    x79 = b + x71
    x80 = 1/x79
    x81 = e*x79
    x82 = x75*(c + x23)
    x83 = x27*x77
    x84 = 1j*x5/4
    x85 = 1j*f
    x86 = a + x85
    x87 = 1/x86
    x88 = f*x86
    x89 = x75*(c + x32)
    x90 = x25*x77
    x91 = b + x85
    x92 = 1/x91
    x93 = f*x91
    x94 = x75*(c + x28)
    x95 = x22*x77
    x96 = 1/(a + x75)
    x97 = 1/(b + x75)
    x98 = 1/(a + x77)
    x99 = 1/(b + x77)
    x100 = -x71
    x101 = a + x100
    x102 = 1/x101
    x103 = -e*x101
    x104 = b + x100
    x105 = 1/x104
    x106 = -e*x104
    x107 = -x85
    x108 = a + x107
    x109 = 1/x108
    x110 = -f*x108
    x111 = b + x107
    x112 = 1/x111
    x113 = -f*x111
    x114 = 1/c
    x115 = 1/x28
    x116 = 1/x32
    x117 = 1/a
    x118 = 1/x24
    x119 = 1/b
    x120 = 1/x21
    x121 = 1/d
    x122 = 1/f
    
    return np.real(x0*(atan(x114*x115*x15) - atan(x114*x116*x3) + atan(x114*x13*x35) - atan(x114*x41*x7)) + x1*atan(x33*x34*x35) + x1*atan(x34*x36*x37) - x1*atan(x34*x38*x39) - x1*atan(x34*x40*x41) - x10*log(b + x28) + x10*log(x0 + x5) - x10*log(x31 + x32) - x11*log(b + x26) + x11*log(x1 + x12) - x11*log(x19 + x31) - x12*(-atan(x118*x121*x3) + atan(x120*x121*x15) + atan(x121*x13*x39) - atan(x121*x37*x7)) - x13*log(x14) + x13*log(x20) + x13*log(x29 + x30) - x15*log(x16) + x15*log(x22) + x15*log(x28 + x29) - x17*log(x2) + x17*log(b + x23) + x17*log(x30 + x31) + x18*log(b + x21) - x18*log(x12 + x5) + x18*log(x24 + x31) + 2*pi*x2 - x3*log(x25) + x3*log(x6) - x3*log(x29 + x32) + x33*log(x42/x43) + x36*log(x44/x45) - x38*log(x46/x47) + x4*(-atan(x10*x116*x117) - atan(x11*x117*x39) + atan(x117*x118*x18) + atan(x117*x17*x35))/2 - x40*log(x48/x49) - x5*(-atan(x115*x122*x40) + atan(x116*x122*x33) - atan(x118*x122*x38) + atan(x120*x122*x36)) - x50*log(x52*(x53 + x54)) + x50*log(x52*(x53 + x55)) + x50*log(x57*(x58 + x59)) - x50*log(x57*(x58 + x60)) + x50*log(x96*(x33 + 1j*x54)) - x50*log(x96*(x33 + 1j*x55)) - x50*log(x97*(x40 + 1j*x59)) + x50*log(x97*(x40 + 1j*x60)) + x61*log(x62*(x63 + x64)) - x61*log(x62*(x63 + x65)) - x61*log(x66*(x67 + x68)) + x61*log(x66*(x67 + x69)) - x61*log(x98*(x38 + 1j*x64)) + x61*log(x98*(x38 + 1j*x65)) + x61*log(x99*(x36 + 1j*x68)) - x61*log(x99*(x36 + 1j*x69)) - x7*log(x27) + x7*log(x9) - x7*log(x23 + x29) - x70*log(-x102*(x103 + x76)) + x70*log(-x102*(x103 + x78)) + x70*log(-x105*(x106 + x82)) - x70*log(-x105*(x106 + x83)) + x70*log(x73*(x74 + x76)) - x70*log(x73*(x74 + x78)) - x70*log(x80*(x81 + x82)) + x70*log(x80*(x81 + x83)) - x8*(-atan(x10*x115*x119) - atan(x11*x119*x37) + atan(x119*x120*x18) + atan(x119*x17*x41))/2 + x84*log(-x109*(x110 + x89)) - x84*log(-x109*(x110 + x90)) - x84*log(-x112*(x113 + x94)) + x84*log(-x112*(x113 + x95)) - x84*log(x87*(x88 + x89)) + x84*log(x87*(x88 + x90)) + x84*log(x92*(x93 + x94)) - x84*log(x92*(x93 + x95)))