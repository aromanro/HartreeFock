// Regression for Eigen allocations crossing the numerical-library boundary.
#include "Basis.h"
#include "RestrictedHartreeFock.h"
#include "UnrestrictedHartreeFock.h"
#include "RestrictedCCSD.h"
#include "RestrictedConfigurationIInteractionSingles.h"
#include <iostream>
#include <memory>
#include <stdexcept>

int main(int argc, char** argv)
{
    try {
        std::cout << std::unitbuf << "Eigen alignment=" << EIGEN_DEFAULT_ALIGN_BYTES
            << ", system allocator=" << EIGEN_MALLOC_ALREADY_ALIGNED << '\n';
        // Client allocates; library destructor frees. The old AVX2/SSE mix
        // crashed here with STATUS_HEAP_CORRUPTION in x64 Release.
        for (int n=1; n<=32; ++n) {
            auto hf=std::make_unique<HartreeFock::RestrictedHartreeFock>();
            hf->h=Eigen::MatrixXd::Ones(n,n);
        }
        Chemistry::Basis basis; basis.Load(argc>1 ? argv[1] : "sto3g.txt");
        Systems::Molecule molecule;
        for (const auto& atom : basis.atoms) if (atom.Z==1) {
            molecule.atoms={atom,atom}; break;
        }
        if (molecule.atoms.size()!=2) throw std::runtime_error("Hydrogen basis not found");
        molecule.atoms[1].SetPosition(Vector3D<double>(0.,0.,1.4));
        molecule.alphaElectrons=molecule.betaElectrons=1;
        molecule.Init();
        auto run=[&](auto& hf) {
            hf.Init(&molecule);
            const double energy=hf.Calculate();
            const double mp2=hf.CalculateMp2Energy();
            if (!hf.converged || std::abs(energy+1.11671432506)>1E-7 || !std::isfinite(mp2))
                throw std::runtime_error("H2 HF/MP2 calculation failed");
            std::cout << "HF=" << energy << ", MP2=" << mp2 << '\n';
        };
        HartreeFock::RestrictedHartreeFock rhf; run(rhf);
        HartreeFock::RestrictedConfigurationIInteractionSingles cis(&rhf);
        if (!cis.Init() || !cis.getSpinOrbitalCISMatrix().allFinite())
            throw std::runtime_error("CIS boundary test failed");
        // Library allocates; client frees the matrix buffer.
        rhf.C.resize(0,0); rhf.h.resize(0,0);
        HartreeFock::UnrestrictedHartreeFock uhf; run(uhf);
        uhf.Cplus.resize(0,0); uhf.Cminus.resize(0,0);
        HartreeFock::RestrictedCCSD cc; run(cc); cc.InitCC();
        const double ccMP2=cc.MP2EnergyFromt4();
        if (!std::isfinite(ccMP2) || std::abs(ccMP2-cc.CalculateMp2Energy())>1E-8)
            throw std::runtime_error("CCSD initial-amplitude boundary test failed");
        std::cout << "Allocation, RHF/UHF, MP2, CIS, and CCSD initialization passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
