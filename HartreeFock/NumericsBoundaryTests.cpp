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
        const double originalMP2=rhf.GetMP2Energy();int amplitudes=0;
        const double visitedMP2=rhf.VisitMp2Amplitudes([&](int,int,int,int,double t){
            if(!std::isfinite(t))throw std::runtime_error("Non-finite exposed MP2 amplitude");++amplitudes;});
        if(!amplitudes || std::abs(visitedMP2-originalMP2)>1e-12)
            throw std::runtime_error("Shared MP2 amplitude visitor changed the RHF energy");
        HartreeFock::RestrictedConfigurationIInteractionSingles cis(&rhf);
        if (!cis.Init() || !cis.getSpinOrbitalCISMatrix().allFinite())
            throw std::runtime_error("CIS boundary test failed");
        // Library allocates; client frees the matrix buffer.
        rhf.C.resize(0,0); rhf.h.resize(0,0);
        HartreeFock::UnrestrictedHartreeFock uhf; run(uhf);
        if (std::abs(uhf.GetMP2Energy()-originalMP2)>1e-9)
            throw std::runtime_error("UMP2 does not recover the RHF MP2 limit");
        uhf.Cplus.resize(0,0); uhf.Cminus.resize(0,0);
        // Independent ordered spin-orbital MP2 sum, using direct AO contraction
        // rather than the cached transform and unique-pair production loops.
        auto checkOpenShell=[&](int z,int na,int nb) {
            Systems::Molecule open;
            for (const auto& atom:basis.atoms) if(atom.Z==z) {open.atoms={atom};break;}
            if(open.atoms.empty()) throw std::runtime_error("Open-shell basis missing");
            if(z==1) {
                const auto hydrogen=open.atoms.front();open.atoms.assign(5,hydrogen);
                for(int k=0;k<5;++k) open.atoms[k].SetPosition(Vector3D<double>(.12*k*k,.17*(k%2),1.6*k));
            }
            open.alphaElectrons=na;open.betaElectrons=nb;open.Init();
            HartreeFock::UnrestrictedHartreeFock u;u.UseDIIS=true;u.Init(&open);u.Calculate();
            if(!u.converged) throw std::runtime_error("Open-shell HF failed");
            const int n=u.numberOfOrbitals;
            std::vector<int> occ,vir;
            auto eps=[&](int p){return p%2?u.eigenvalsplus(p/2):u.eigenvalsminus(p/2);};
            for(int p=0;p<2*n;++p) ((p/2<(p%2?na:nb))?occ:vir).push_back(p);
            auto coulomb=[&](int p,int a,int q,int b) {
                if(p%2!=a%2 || q%2!=b%2) return 0.;
                const auto& c=p%2?u.Cplus:u.Cminus;
                const auto& d=q%2?u.Cplus:u.Cminus;
                double g=0.;
                for(int mu=0;mu<n;++mu) for(int nu=0;nu<n;++nu)
                    for(int la=0;la<n;++la) for(int si=0;si<n;++si)
                        g+=c(mu,p/2)*c(nu,a/2)*d(la,q/2)*d(si,b/2)
                            *u.integralsRepository.getElectronElectron(mu,nu,la,si);
                return g;
            };
            double reference=0.;
            for(int i:occ) for(int j:occ) for(int a:vir) for(int b:vir) {
                const double g=coulomb(i,a,j,b)-coulomb(i,b,j,a);
                reference+=.25*g*g/(eps(i)+eps(j)-eps(a)-eps(b));
            }
            const double actual=u.CalculateMp2Energy();
            if(std::abs(actual-reference)>1e-10) throw std::runtime_error("UMP2 differs from direct AO spin-orbital reference");
            double emitted[3]={};
            u.VisitMp2Amplitudes([&](bool up,bool tau,int i,int j,int a,int b,double t){
                const double den=(up?u.eigenvalsplus:u.eigenvalsminus)(i)
                    +(tau?u.eigenvalsplus:u.eigenvalsminus)(j)
                    -(up?u.eigenvalsplus:u.eigenvalsminus)(a)
                    -(tau?u.eigenvalsplus:u.eigenvalsminus)(b);
                emitted[up==tau?(up?0:1):2]+=t*t*den;
            });
            if(std::abs(emitted[0]+emitted[1]+emitted[2]-reference)>1e-10)
                throw std::runtime_error("UMP2 amplitude/energy inconsistency");
            if(z==1 && !(emitted[0]<-1e-8 && emitted[1]<-1e-8 && emitted[2]<-1e-8))
                throw std::runtime_error("H5 fixture must exercise all UMP2 spin channels");
            std::swap(u.Cplus,u.Cminus);std::swap(u.eigenvalsplus,u.eigenvalsminus);
            std::swap(u.occupiedPlus,u.occupiedMinus);
            if(std::abs(u.CalculateMp2Energy()-reference)>1e-10)
                throw std::runtime_error("UMP2 spin-swap invariance failed");
            std::cout<<"UMP2 Z="<<z<<" direct="<<reference<<" aa="<<emitted[0]
                <<" bb="<<emitted[1]<<" ab="<<emitted[2]<<'\n';
        };
        checkOpenShell(3,2,1);checkOpenShell(5,3,2);checkOpenShell(1,3,2);
        molecule.alphaElectrons=1;molecule.betaElectrons=0;molecule.Init();
        HartreeFock::UnrestrictedHartreeFock one;one.Init(&molecule);one.Calculate();
        int emitted=0;
        if(!one.converged || one.VisitMp2Amplitudes([&](bool,bool,int,int,int,int,double){++emitted;})!=0. || emitted)
            throw std::runtime_error("One-electron UMP2 must vanish");
        molecule.alphaElectrons=molecule.betaElectrons=1;molecule.Init();
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
