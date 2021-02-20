#!/bin/bash
#these diff commands should yeild no difference

diff geopack_08/try_using_1_subroutine/ZZZrun.sh mead/try_using_1_subroutine/ZZZrun.sh
diff geopack_08/try_using_2_subroutine/ZZZrun.sh mead/try_using_2_subroutine/ZZZrun.sh

diff mead/try_using_1_subroutine/main.py mead/try_using_2_subroutine/main.py
diff geopack_08/try_using_1_subroutine/main.py geopack_08/try_using_2_subroutine/main.py
