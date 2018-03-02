# Copyright 1999-2018 Gentoo Foundation
# Distributed under the terms of the GNU General Public License v2

EAPI=6

inherit cmake-utils

DESCRIPTION="libscientific for linear algebra calculations"
HOMEPAGE="https://github.com/gmrandazzo/libscientific"


if [[ ${PV} == 9999 ]]; then
	inherit git-r3
	EGIT_REPO_URI="https://github.com/gmrandazzo/${PN}.git"
else
    SRC_URI="https://github.com/gmrandazzo/${PN}/archive/v1.0.0.tar.gz -> ${P}.tar.gz"
	KEYWORDS="amd64 x86"
fi

LICENSE="LGPLv3"
SLOT="0"
KEYWORDS="~amd64 ~x86"
IUSE=""

DEPEND="dev-util/cmake"
RDEPEND="${DEPEND}"

src_install() {
	cmake-utils_src_install
}

