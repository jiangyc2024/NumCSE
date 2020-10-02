

# Tests of autofocus.hpp

> There are five functions in autofocus.hpp:
* `save_image(double focus)`
* `void plot_freq(double focus)`
* `high_frequency_content(const MatrixXd & M)`
* `plotV()`
* `autofocus()`

> The third one `high_frequency_content` is only used for plotV(). There is no test for this.

> The first two functions are defined on the same inputs, namely `focus = 0, 1, 2, 3`.