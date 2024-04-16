# pathIntegral (from Diwina project)
# based on FeeLLGood – A micromagnetic solver ![Build Status](https://github.com/feellgood/FeeLLGood/actions/workflows/tests.yml/badge.svg)

pathIntegral is a STXM image simulator using feeLLGood simulations inputs. feeLLGood is a micromagnetic solver using finite element technique to integrate Landau Lifshitz Gilbert equation, developped by JC Toussaint & al. The code is being modified without any warranty it works. A dedicated website can be found [here][]  

### Dependencies

* C++17 and the STL
* [GMSH][]
* [TBB][]
* [libpng][]
* [yaml-cpp][]
* [Duktape][] 2.7.0
* [Eigen][] ≥ 3.3

### License

Copyright (C) 2024  Jean-Christophe Toussaint, C. Thirion and E. Bonet.

pathIntegral is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

Additional permission under GNU GPL version 3 section 7: If you modify this Program, or any covered work, by linking or combining it with the Intel® MKL library (or a modified version of that library), containing parts covered by the terms of Intel Simplified Software License, the licensors of this Program grant you additional permission to convey the resulting work.

The libraries used by pathIntegral are distributed under different licenses, and this is documented in their respective Web sites.

[here]: https://feellgood.neel.cnrs.fr/
[GMSH]: https://gmsh.info/
[TBB]: https://www.threadingbuildingblocks.org/
[libpng]: http://libpng.org/pub/png/libpng.html
[yaml-cpp]: https://github.com/jbeder/yaml-cpp
[Duktape]: https://duktape.org/
[Eigen]: https://eigen.tuxfamily.org/
