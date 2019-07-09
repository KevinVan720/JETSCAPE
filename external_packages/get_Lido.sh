#!/usr/bin/env bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
#
# Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
#
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

# download the code package
git clone -b Lido-2.0 https://github.com/keweiyao/Duke-Lido.git Lido
cd Lido
git checkout a3036d1c3ccebf6631c4d61615c98100c17f7938
cd ..


