#pragma once

#include <Eigen/Dense>

/* \brief set_focus Given a double, retuns an image "taken" with
 * the double as "focusing parameter".
 * The focus parameter "f0" simulate the focal length chosen by
 * a digital camera. The function returns a $n \times m$ Matrix of doubles,
 * whose
 * values represent a grayscale image $n \times m$, with values
 * between 0 and 255, and where 0 represents black, and 255 represents
 * white.
 * \param[in] f0 The focal length parameter.
 * \return A MatrixXd, which contains the grey scale values of the image.
 */
Eigen::MatrixXd set_focus(double f0);

#include "._SUPER_SECRET_FILE.hpp"
