def bounded_first_derivative(xs, ys, apply_bounds=True):
    '''
    intput 3x points
    output 1st-order deriv appx
    '''
    dxp = xs[2] - xs[1]
    dp = (ys[2] - ys[1]) / dxp
    dxm = xs[1] - xs[0]
    dm = (ys[1] - ys[0]) / dxm
    if dp * dm < 0:
        return 0.0
    di = (dxm * dp + dxp * dm) / (dxp + dxm)
    dimin = min(abs(dp), abs(dm))
    if apply_bounds:
        return max(-3*dimin, min(3*dimin, di))
    else:
        return di

def monotone_cubic_interpolation(xs, ys, xt, apply_bounds=True):
    '''
    assume len(xs) >= 3
    xs is monotonically increasing
    ys is monotonic
    '''
    EPS = 1e-6
    def phi(t):
        return 3*t*t - 2*t*t*t
    def psi(t):
        return t*t*t - t*t
    # locate i s.t. xt inside (x[i], x[i+1])
    if xt < xs[0] + EPS:
        return ys[0]
    if xt > xs[-1] - EPS:
        return ys[-1]
    n = len(xs)
    for i in range(n):
        if xt < xs[i]:
            i = i - 1
            break
    # 1st-order derivatives
    # d_i
    if i == 0:
        di = (ys[1] - ys[0]) / (xs[1] - xs[0])
    else:
        di = bounded_first_derivative(xs[i-1:], ys[i-1:], apply_bounds)
    # d_{i+1}
    if i == n-2:
        dip1 = (ys[n-1] - ys[n-2]) / (xs[n-1] - xs[n-2])
    else:
        dip1 = bounded_first_derivative(xs[i:], ys[i:], apply_bounds)
    # formula (2) in reference 1
    dx = xs[i+1] - xs[i]
    dxp = (xs[i+1] - xt) / dx
    dxm = 1 - dxp
    return (ys[i] * phi(dxp) + ys[i+1] * phi(dxm) + di * (-dx) * psi(dxp) + dip1 * dx * psi(dxm))
    
    # wiki page example

xs = [0.1*(i+1) for i in range(9)]
ys = [3, 2.8, 2.5, 1, 0.95, 0.7, 0.5, 0.2, 0.1]

mxs = [0.01*(i+10) for i in range(81)]
uys = [monotone_cubic_interpolation(xs, ys, xt, False) for xt in mxs]
mys = [monotone_cubic_interpolation(xs, ys, xt) for xt in mxs]

import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('ggplot')
mpl.rcParams['figure.figsize'] = [12, 6]

fig, ax = plt.subplots()

ax.plot(xs, ys, marker='s', linewidth=0, label='Data points')
ax.plot(mxs, uys, label='Cubic')
ax.plot(mxs, mys, label='Monotone cubic')
ax.legend(facecolor='white')
plt.show()
