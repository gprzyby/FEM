from math import sqrt

#define 2D shape functions
#eta = η
#ksi = ξ
shape1d_fun = (lambda ksi:  0.5 * (1 - ksi), \
                 lambda ksi:  0.5 * (ksi + 1))

shape1d_dksi = (-0.5, \
                  0.5)

shape2d_fun = (lambda ksi, eta: 0.25 * (1 - ksi)*(1 - eta), \
                 lambda ksi, eta: 0.25 * (1 + ksi)*(1 - eta), \
                 lambda ksi, eta: 0.25 * (1 + ksi)*(1 + eta), \
                 lambda ksi, eta: 0.25 * (1 - ksi)*(1 + eta))

shape2d_dksi = (lambda ksi, eta: -0.25 * (1 - eta), \
                  lambda ksi, eta: 0.25 * (1 - eta), \
                  lambda ksi, eta: 0.25 * (1 + eta), \
                  lambda ksi, eta: -0.25 * (1 + eta))

shape2d_deta = (lambda ksi, eta: -0.25 * (1 - ksi), \
                  lambda ksi, eta: -0.25 * (1 + ksi), \
                  lambda ksi, eta: 0.25 * (1 + ksi), \
                  lambda ksi, eta: 0.25 * (1 - ksi))

#define points for gauss integral
gauss_points = (
                [
                    (0., 2.)
                ],
                (
                    (-1./(3.**0.5), 1.),
                    (1./(3.**0.5), 1.)
                ),
                (
                    (-(3./5.)**0.5, 5./9.),
                    (0, 8./9.),
                    ((3./5.)**0.5, 5./9.),
                ),
                (
                    (-sqrt(3./7. + (2./7.)*sqrt(6./5.)), (18. - sqrt(30.))/36.),
                    (-sqrt(3./7. - (2./7.)*sqrt(6./5.)), (18. + sqrt(30.))/36.),
                    (sqrt(3./7. - (2./7.)*sqrt(6./5.)), (18. + sqrt(30.))/36.),
                    (sqrt(3./7. + (2./7.)*sqrt(6./5.)), (18. - sqrt(30.))/36.)
                ),

)

