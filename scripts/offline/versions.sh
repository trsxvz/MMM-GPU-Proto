# Versions packed by scripts/offline/bundle.sh: the ones of the spack
# environments of the prototype (deps/spack/README.md) and the base of the
# patches of MGIS (scripts/sync-mgis.sh).
PROTO_SPACK_VERSION=v1.2.2
PROTO_SPACK_PACKAGES_COMMIT=3b7fdffc24
PROTO_MGIS_URL=https://github.com/thelfer/MFrontGenericInterfaceSupport.git
PROTO_MGIS_BASE=b1afeb3
# MFEM built from its sources by `deploy.sh mfem` on the machines whose
# libraries come from their modules (no spack environment, see
# docs/install.md): archive of the release (URL and checksum of the recipe
# of spack)
PROTO_MFEM_VERSION=4.10.0
PROTO_MFEM_URL=https://bit.ly/4iEAcbE
PROTO_MFEM_SHA256=d01c26662ab96042ec5e3443a58ee9eb3c5c78c6632d899240b6e485f3a414cf
