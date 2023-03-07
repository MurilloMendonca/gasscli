# Maintainer: Murillo Ventura <murillo.ventura1711@gmail.com>

pkgname=gasscli
pkgver=1.0
pkgrel=1
pkgdesc="A CLI tool for the newGASS genetic algorithm library"
arch=('x86_64')
url="https://github.com/MurilloMendonca/gasscli"
license=('MIT')
depends=('boost' 'curl')
source=("git+$url")
md5sums=('SKIP')

build() {
  cd "$srcdir/gasscli"
  git submodule update --init --recursive
  make
}

package() {
  cd "$srcdir/gasscli"
  install -Dm755 gasscli "$pkgdir"/usr/bin/gasscli
  # Call configure.sh to create the cache folder
  ./configure.sh
}

