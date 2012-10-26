//-------------------------------------------------------------------
// Copy this file into your HEPMCDIR/test/ and compile it
//-------------------------------------------------------------------
//

#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#ifndef HEPMC_IO_ASCII_REMOVED
#include "HepMC/IO_Ascii.h"
#endif
#ifdef HEPMC_HAS_IO_GENEVENT
#include "HepMC/IO_GenEvent.h"
#endif
#include "HepMC/IO_AsciiParticles.h"

// define methods and classes used by this test
#include "IsGoodEvent.h"
#include "testHepMCMethods.h"
#include <sstream>
#include <list>
#include <cmath>


void find_protons();
void my_main();
std::vector<HepMC::GenParticle*> get_in_particles(HepMC::GenEvent* evt);
std::vector<HepMC::GenParticle*> get_interm_particles(HepMC::GenEvent* evt);
std::vector<HepMC::GenParticle*> get_out_particles(HepMC::GenEvent* evt);
void print_latex(HepMC::GenEvent* evt);
void print_latex_in_particles(std::vector <HepMC::GenParticle*> in_p);
void print_latex_beam_particles(std::vector <HepMC::GenParticle*> beam_p);
void print_latex_interm_particles(std::vector <HepMC::GenParticle*> interm_p,std::map<HepMC::GenVertex *, int> level_map, HepMC::GenEvent * evt);
void print_latex_out_particles(std::vector <HepMC::GenParticle*> out_p, std::map<HepMC::GenVertex *, int> level_map, HepMC::GenEvent * evt);
std::string  pdg_to_line(int pdg_id);
HepMC::GenEvent  strip_event(HepMC::GenEvent* evt);
void recur_copy_particles(HepMC::GenVertex * v_in, HepMC::GenEvent* evt);
void add_protons(HepMC::GenEvent *evt);
void add_interaction(HepMC::GenEvent *evt, HepMC::GenEvent *evt_in );
HepMC::GenVertex * find_interaction_vertex(HepMC::GenEvent * evt);
std::string pdg_to_label(int pdg_id);
bool particle_has_ancestor(HepMC::GenParticle * p, HepMC::GenParticle * comp_p);
bool interesting_pdg_id(int pdg_id);
bool interesting_end_vertex(HepMC::GenParticle * p);
void recursive_sort_out_particles(HepMC::GenVertex * curr_v, int curr_lev, 
                                        std::map<HepMC::GenVertex *, int > & v_map, std::vector<HepMC::GenParticle* > & p_out );

int main() { 
    my_main();
    return 0;
}

void my_main(){
    HepMC::IO_GenEvent my_events( "./pythia-output.hepmc" ,std::ios::in);
//    HepMC::IO_GenEvent my_events( "./testFlow.out3" ,std::ios::in);
//    HepMC::IO_GenEvent my_events( "./testIOGenEvent.input" ,std::ios::in);
    my_events.use_input_units( HepMC::Units::GEV, HepMC::Units::MM );
    HepMC::GenEvent* evt = my_events.read_next_event();
////////////////// WHILE lOOP 
    int count =0;
    while ( evt ) {
//        std::cout << std::endl << std::endl << "Event: " << count << std::endl << std::endl;
        count ++;
///////////////////////////
        if (count==25){
            evt->print();
            HepMC::GenEvent* strip_evt = new HepMC::GenEvent(strip_event(evt));
            print_latex(strip_evt);
        }
//////////////// WHILE LOOP ///////////////
	    delete evt;
	    my_events >> evt;
    }
//////////////// WHILE LOOP ///////////////
}




bool interesting_pdg_id(int pdg_id){
    std::map<int,int> yes;
    yes[ 5]=0;
    yes[ 6]=0;
    yes[15]=0;
    yes[23]=0;
    yes[24]=0;
    yes[25]=0;
    yes[35]=0;
    yes[36]=0;
    yes[37]=0;
    if (abs(pdg_id) > 1000000) return true;
    else if (yes.count(abs(pdg_id))) return true;
    else return false;
}

bool interesting_end_vertex(HepMC::GenParticle * p){
    // check whether particle is supersymmetric (or 3rd generation, Higgs,  W or Z ) and has a decay vertex
    return (interesting_pdg_id(p->pdg_id()) && p->end_vertex() );
}

void recur_copy_particles(HepMC::GenVertex * v_in, HepMC::GenVertex * v_curr , HepMC::GenEvent* evt){

    for (HepMC::GenVertex::particles_out_const_iterator p_in = v_in->particles_out_const_begin(); 
         p_in != v_in->particles_out_const_end(); p_in++  ){
//      create new HepMC::GenParticle
        HepMC::GenParticle * part = new HepMC::GenParticle();
//      VERY IMPORTANT: shallow copy the content of **p_in to *part
        *part = **p_in;
//      add this to the current vertex
        v_curr->add_particle_out(part);

        if (interesting_end_vertex((*p_in))){
            HepMC::GenVertex * v = new HepMC::GenVertex();
            evt->add_vertex(v);
            v->add_particle_in(part);
            recur_copy_particles((*p_in)->end_vertex(), v ,evt);
        }
    }
}

void add_protons(HepMC::GenEvent *evt){
    //  CREATE these particles and vertices:
    //   p1 - pv1 - p1o
    //         
    //         
    //         
    //   p2 - pv2 - p2o
    // initialise Vertices and particles
    HepMC::GenVertex * pv1= new HepMC::GenVertex();
    HepMC::GenVertex * pv2= new HepMC::GenVertex();
    HepMC::GenParticle* p1 = new HepMC::GenParticle( HepMC::FourVector(0,0,7000,7000), 2212);
    HepMC::GenParticle* p2 = new HepMC::GenParticle( HepMC::FourVector(0,0,-7000,7000), 2212);
    HepMC::GenParticle* p1o = new HepMC::GenParticle( HepMC::FourVector(0,0,7000,7000), 2212);
    HepMC::GenParticle* p2o = new HepMC::GenParticle( HepMC::FourVector(0,0,-7000,7000), 2212);
    // add particles to vertices
    pv1->add_particle_in(p1);
    pv1->add_particle_out(p1o);
    pv2->add_particle_in(p2);
    pv2->add_particle_out(p2o);
    // add vertices to event
    evt->add_vertex(pv1);
    evt->add_vertex(pv2);
    // set beam particles 
    evt->set_beam_particles(p1,p2);
}

bool particle_has_ancestor(HepMC::GenParticle * p, HepMC::GenParticle * comp_p){
    HepMC::GenVertex * v = p->production_vertex();
    bool match=false;
	for ( HepMC::GenVertex::particle_iterator part = v->particles_begin(HepMC::ancestors);
	      part != v->particles_end(HepMC::ancestors); 
	      ++part ) {
        if( (*part)==comp_p){ match=true;break; }
    }
    return match;
}

HepMC::GenVertex * find_interaction_vertex(HepMC::GenEvent * evt){
    //create a vector of potential interaction vertices, i.e. with 2 incoming particles
    std::vector<HepMC::GenVertex*> pot_iv;
    std::vector<HepMC::GenVertex*> cand_iv;

    for ( HepMC::GenEvent::vertex_iterator v = evt->vertices_begin(); v != evt->vertices_end(); ++v ) {
        if((*v)->particles_in_size()==2){ 
            pot_iv.push_back((*v));
        }
    }
    // loop over potential vertices 
    for (unsigned int i = 0; i < pot_iv.size(); i++){
        // check whether each of the potientially interacting particles came from a different beam particle
        HepMC::GenParticle * p1 = evt->beam_particles().first;
        HepMC::GenParticle * p2 = evt->beam_particles().second;
        HepMC::GenParticle *ip1 = *(pot_iv[i]->particles_in_const_begin()); 
        HepMC::GenParticle *ip2 = *(++(pot_iv[i]->particles_in_const_begin())  ); 
        if(ip1->production_vertex() && ip2->production_vertex()){
            ip1->print();
            ip2->print();
            bool sp1_from_p1= particle_has_ancestor( ip1,p1 ); 
            bool sp1_from_p2= particle_has_ancestor( ip1,p2 ); 
            bool sp2_from_p1= particle_has_ancestor( ip2,p1 ); 
            bool sp2_from_p2= particle_has_ancestor( ip2,p2 ); 
            if ( ((sp1_from_p1 && sp2_from_p2  ) || (sp1_from_p2 && sp2_from_p1)   )  ) cand_iv.push_back( pot_iv[i]) ;
        }
    }
    if(cand_iv.size()==1){ 
        std::cout << "Interaction vertex found:" << std::endl;
        cand_iv[0]->print();
        return cand_iv[0];
    }
    else if(cand_iv.size()>1) {
        std::cerr << "Ambigous interaction vertices found:" << std::endl;
        for (unsigned int i = 0 ; i< cand_iv.size(); i++){
            std::cerr << "Potential vertex candidate " << i << std::endl;
            cand_iv[i]->print(std::cerr);
        }
        return NULL;
    }
    else {
        std::cerr << "No interaction vertex found" << std::endl;
        return NULL;
    }
}

void add_interaction(HepMC::GenEvent *evt, HepMC::GenEvent *evt_in ){
    //  CREATE vi and connect to pv1 and pv2
    //  vi is from the interaction vertex of evt_in
    //  the two incoming particles are also from evt_in
    //   p1 - pv1 - p1o
    //         \
    //         vi
    //         /
    //   p2 - pv2 - p2o
    HepMC::GenVertex * iv_in = evt_in->signal_process_vertex();
    if (!iv_in) iv_in=find_interaction_vertex(evt_in);

    HepMC::GenParticle * p1 = new HepMC::GenParticle();
    HepMC::GenParticle * p2 = new HepMC::GenParticle();

    //get the incoming particles from the interaction vertex
    *p1 = *( *(iv_in->particles_in_const_begin()) );
    *p2 = *( *( ++(iv_in->particles_in_const_begin()) ) );

    HepMC::GenVertex * iv = new HepMC::GenVertex();
    iv->add_particle_in(p1);
    iv->add_particle_in(p2);
    if (evt->valid_beam_particles()){
        (evt->beam_particles().first)->end_vertex()->add_particle_out(p1);
        (evt->beam_particles().second)->end_vertex()->add_particle_out(p2);
    }
    evt->add_vertex(iv );
    evt->set_signal_process_vertex(iv);
}

HepMC::GenEvent  strip_event(HepMC::GenEvent* in_evt){
    // create new event that will contain the stripped version of in_evt
    HepMC::GenEvent* evt = new HepMC::GenEvent( 0,0  );
    evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);
     
    // add the incomming protons + vertices and the interaction vertex
    add_protons(evt);
    add_interaction(evt,in_evt);

    // add everything that follows the interaction vertex
    recur_copy_particles(in_evt->barcode_to_vertex(-5),evt->signal_process_vertex(),evt);    

    return (* evt);
}

std::vector<HepMC::GenParticle*> get_in_particles(HepMC::GenEvent* evt){
        std::vector<HepMC::GenParticle*> in_p;
        for (HepMC::GenEvent::particle_iterator p = evt->particles_begin(); 
             p != evt->particles_end(); p++  ){
            if(!(*p)->production_vertex() && (*p)->end_vertex()){
                in_p.push_back((*p));
            }
        }
        
        return in_p;
}

std::vector<HepMC::GenParticle*> get_beam_particles(HepMC::GenEvent* evt){
        std::vector<HepMC::GenParticle*> beam_p;
        if (evt->valid_beam_particles()){
            beam_p.push_back(evt->beam_particles().first);
            beam_p.push_back(evt->beam_particles().second);
        }
        else{
            std::cout << "WARNING: no valid beam particles" << std::endl;
            std::cerr << "WARNING: no valid beam particles" << std::endl;
        }
        return beam_p;
}

std::vector<HepMC::GenParticle*> get_interm_particles(HepMC::GenEvent* evt){
        std::vector<HepMC::GenParticle*> interm_p;
        for (HepMC::GenEvent::particle_iterator p = evt->particles_begin(); 
             p != evt->particles_end(); p++  ){
            if((*p)->production_vertex() && (*p)->end_vertex() ){
                interm_p.push_back((*p));
            }
        }
        
        return interm_p;
}

std::vector<HepMC::GenParticle*> get_out_particles(HepMC::GenEvent* evt){
        std::vector<HepMC::GenParticle*> out_p;
        for (HepMC::GenEvent::particle_iterator p = evt->particles_begin(); 
             p != evt->particles_end(); p++  ){
            if(!(*p)->end_vertex()){
                out_p.push_back((*p));
            }
        }
        
        return out_p;
}

std::string  pdg_to_line(int pdg_id){
    if (pdg_id<0) pdg_id=-pdg_id;
    std::string s;
    // lepton or quark
    if ((pdg_id >0 && pdg_id <9 ) || (pdg_id >10 && pdg_id <18)) s = "fermion";
    // gluon
    else if (pdg_id == 21) s="gluon";
    // gamma, W, Z
    else if (pdg_id >22 && pdg_id < 25) s="boson";
    // proton
    else if (pdg_id==2212) s="dbl_plain_arrow";
    // neutralinos and charginos
    else if (pdg_id >1000021 && pdg_id < 1000037) s="gaugino";
    // gluino
    else if (pdg_id == 1000021) s = "gluino";
    else s="dashes";
    return s;
}


std::string pdg_to_label(int pdg_id){
    std::map<int , std::string> pdgmap;
    // quarks and anti-quarks
    pdgmap[1]="d";
    pdgmap[2]="u";
    pdgmap[3]="s";
    pdgmap[4]="c";
    pdgmap[5]="b";
    pdgmap[6]="t";
    pdgmap[-1]="\\bar{d}";
    pdgmap[-2]="\\bar{u}";
    pdgmap[-3]="\\bar{s}";
    pdgmap[-4]="\\bar{c}";
    pdgmap[-5]="\\bar{b}";
    pdgmap[-6]="\\bar{t}";
    // leptons
    pdgmap[11]="e^-";
    pdgmap[13]="\\mu^-";
    pdgmap[15]="\\tau^-";
    pdgmap[-11]="e^+";
    pdgmap[-13]="\\mu^+";
    pdgmap[-15]="\\tau^+";
    pdgmap[12]="\\nu_{e}";
    pdgmap[14]="\\nu_{\\mu }";
    pdgmap[16]="\\nu_{\\tau }";
    pdgmap[-12]="\\bar{\\nu}_{e}";
    pdgmap[-14]="\\bar{\\nu}_{\\mu }";
    pdgmap[-16]="\\bar{\\nu}_{\\tau }";
    // gauge bosons
    pdgmap[21]="g";
    pdgmap[22]="\\gamma";
    pdgmap[23]="Z";
    pdgmap[24]="W^+";
    pdgmap[-24]="W^-";
    pdgmap[25]="h^0";
    pdgmap[35]="H^0";
    pdgmap[36]="A^0";
    pdgmap[37]="H^+";
    pdgmap[-37]="H^-";
    // proton
    pdgmap[2212]="p";
    pdgmap[-2212]="\\bar{p}";
    // SUSY particles
    // quarks 
    pdgmap[1000001]="\\tilde{d}_L";
    pdgmap[1000002]="\\tilde{u}_L";
    pdgmap[1000003]="\\tilde{s}_L";
    pdgmap[1000004]="\\tilde{c}_L";
    pdgmap[1000005]="\\tilde{b}_1";
    pdgmap[1000006]="\\tilde{t}_1";
    pdgmap[2000001]="\\tilde{d}_R";
    pdgmap[2000002]="\\tilde{u}_R";
    pdgmap[2000003]="\\tilde{s}_R";
    pdgmap[2000004]="\\tilde{c}_R";
    pdgmap[2000005]="\\tilde{b}_2";
    pdgmap[2000006]="\\tilde{t}_2";
    pdgmap[-1000001]="\\bar{\\tilde{d}}_L";
    pdgmap[-1000002]="\\bar{\\tilde{u}}_L";
    pdgmap[-1000003]="\\bar{\\tilde{s}}_L";
    pdgmap[-1000004]="\\bar{\\tilde{c}}_L";
    pdgmap[-1000005]="\\bar{\\tilde{b}}_1";
    pdgmap[-1000006]="\\bar{\\tilde{t}}_1";
    pdgmap[-2000001]="\\bar{\\tilde{d}}_R";
    pdgmap[-2000002]="\\bar{\\tilde{u}}_R";
    pdgmap[-2000003]="\\bar{\\tilde{s}}_R";
    pdgmap[-2000004]="\\bar{\\tilde{c}}_R";
    pdgmap[-2000005]="\\bar{\\tilde{b}}_2";
    pdgmap[-2000006]="\\bar{\\tilde{t}}_2";
    // sleptons
    pdgmap[1000011]="\\tilde{e}_L^-";
    pdgmap[1000013]="\\tilde{\\mu}_L^-";
    pdgmap[1000015]="\\tilde{\\tau}_1^-";
    pdgmap[2000011]="\\tilde{e}_R^-";
    pdgmap[2000013]="\\tilde{\\mu}_R^-";
    pdgmap[2000015]="\\tilde{\\tau}_2^-";
    pdgmap[1000012]="\\tilde{\\nu}_{eL}";
    pdgmap[1000014]="\\tilde{\\nu}_{\\mu eL}";
    pdgmap[1000016]="\\tilde{\\nu}_{\\tau eL}";
    pdgmap[-1000011]="\\tilde{e}_L^+";
    pdgmap[-1000013]="\\tilde{\\mu}_L^+";
    pdgmap[-1000015]="\\tilde{\\tau}_1^+";
    pdgmap[-2000011]="\\tilde{e}_R^+";
    pdgmap[-2000013]="\\tilde{\\mu}_R^+";
    pdgmap[-2000015]="\\tilde{\\tau}_2^+";
    pdgmap[-1000012]="\\bar{\\tilde{\\nu}}_{eL}";
    pdgmap[-1000014]="\\bar{\\tilde{\\nu}}_{\\mu L}";
    pdgmap[-1000016]="\\bar{\\tilde{\\nu}}_{\\tau L}";

    // gluino
    pdgmap[1000021]="\\tilde{g}";
    // neutralinos
    pdgmap[1000022]="\\tilde{\\chi}^0_1";
    pdgmap[1000023]="\\tilde{\\chi}^0_2";
    pdgmap[1000025]="\\tilde{\\chi}^0_3";
    pdgmap[1000026]="\\tilde{\\chi}^0_4";
    // charginos
    pdgmap[1000024]="\\tilde{\\chi}^+_1";
    pdgmap[-1000024]="\\tilde{\\chi}^-_1";
    pdgmap[1000037]="\\tilde{\\chi}^+_2";
    pdgmap[-1000037]="\\tilde{\\chi}^-_2";

    if (pdgmap.count(pdg_id)) return pdgmap[pdg_id];
    else{ 
        std::stringstream pdg_str;
        pdg_str << pdg_id;
        return ("pdg"+pdg_str.str());
    }
}

void recursive_sort_out_particles(HepMC::GenVertex * curr_v, int curr_lev, 
                                        std::map<HepMC::GenVertex *, int > & v_map, std::vector<HepMC::GenParticle* > & p_out ){
    // The idea is to loop recursively like in recur_copy_particles(), keeping track of
    //  * the level of the vertex, i.e. how far in the decay chain
    //  * in which order the out_particles are reached. This is also the order in which they are printed
    v_map[curr_v]=curr_lev;
    for (HepMC::GenVertex::particle_iterator it = curr_v->particles_begin(HepMC::children); 
                it != curr_v->particles_end(HepMC::children); it++){
        if ((*it)->end_vertex()) recursive_sort_out_particles((*it)->end_vertex(),curr_lev++,v_map, p_out);
        else p_out.push_back(*it);
    }
}

std::pair< std::vector<HepMC::GenParticle*> , std::map<HepMC::GenVertex*, int> > get_sorted_out_particles(  HepMC::GenEvent* evt){
    // Event should look like
    //             p1o
    //            /   
    //   p1 - pv1         smp1
    //         \    sp1<  sp1  
    //         iv<
    //         /    sp2<  sp2
    //   p2 - pv2         smp2  
    //            \
    //             p2o
    // the outgoing vertices should come in this order: p1o, smp1, sp1,  sp2, smp2 , p2o
    // Note: p1o has beam particle p1 (and not p2) in its ancestors, p2o has p2 (and not p1), whereas the rest have both.
    // The rest we can sort recursively


    // get unsorted our particles
    std::vector<HepMC::GenParticle*> p_in=get_out_particles(evt);

    // initiate level map
    std::map<HepMC::GenVertex *, int >  level_map;

    // get the interaction vertex
    HepMC::GenVertex * iv = evt->signal_process_vertex();
    if (!iv) iv=find_interaction_vertex(evt);
    if (!iv) std::cerr << "WARNING: This event has no interaction vertex" << std::endl;

    // only start the whole chain if the event has well defined beam particles
    if (evt->valid_beam_particles()){
        std::vector<HepMC::GenParticle*> p_sort;
        // define 3 categories pointing back to: 1)p1  2) interaction vertec 3) p2 
        std::vector<HepMC::GenParticle*> p1_l;
        std::vector<HepMC::GenParticle*> p2_l;
        std::vector<HepMC::GenParticle*> rest;

        // select particles p1o and p2o
        for (unsigned int i=0 ; i<p_in.size() ; i++){
            int p1=0;
            int p2=0;
            if( particle_has_ancestor( p_in[i] ,evt->beam_particles().first )){p1++;} 
            if( particle_has_ancestor( p_in[i] ,evt->beam_particles().second)){p2++;}

            if( p1==1 && p2==0){p1_l.push_back(p_in[i]); }
            else if( p1==0 && p2==1){p2_l.push_back(p_in[i]); }
        }

        // get the rest of the out particles that originated from the interaction vertex, completely sorted, 
        // meanwhile, keep track of the level of the vertices 
    //    std::map<HepMC::GenVertex *, int >  level_map;
        recursive_sort_out_particles(iv,0,level_map,rest );

        // fill vector in correct order
        for(std::vector<HepMC::GenParticle*>::iterator it=p1_l.begin(); it!= p1_l.end(); it++){p_sort.push_back(*it ); }
        for(std::vector<HepMC::GenParticle*>::iterator it=rest.begin(); it!= rest.end(); it++){p_sort.push_back(*it ); }
        for(std::vector<HepMC::GenParticle*>::iterator it=p2_l.begin(); it!= p2_l.end(); it++){p_sort.push_back(*it ); }
        return std::make_pair (p_sort,level_map);
    }
    else {
            // make all vertices descendent from the interaction vertex level 0
            for (HepMC::GenVertex::vertex_iterator it = iv->vertices_begin(HepMC::descendants); 
                it != iv->vertices_end(HepMC::descendants); it++){
                level_map[(*it)]=0;
            }
            return std::make_pair(p_in,level_map);
    }
}

void print_latex_in_particles(std::vector <HepMC::GenParticle*> in_p){
        // label the incoming particles "iv"+str(p.barcode), e.g. iv1 if the barcode is 1
        // output will look like '\fmfleft{iv1,iv3,}'
        std::stringstream in_vertices;
        in_vertices << "\\fmfleft{" ;
        for (unsigned int i=0; i< in_p.size(); i++){
            in_vertices << "iv" << in_p[i]->barcode(); 
            if (i< (in_p.size()-1) ) in_vertices << ","; 
        }
        in_vertices << "}" ;
        std::cout << in_vertices.str() << std::endl;
        // draw the lines: begin vertices are given by code above, the inner vertices are given by
        // "v"+endvertex(barcode)
        for (unsigned int i = 0; i< in_p.size(); i++){
            std::cout << "\\fmf{" 
                      << pdg_to_line(in_p[i]->pdg_id()) 
                      << ",lab=$" << pdg_to_label(in_p[i]->pdg_id()) << "$"  
                      <<"}{"
                      << "iv" << in_p[i]->barcode() << ","
                      << "v" << -(in_p[i]->end_vertex()->barcode())
                      << "}" << std::endl;
        }

}

void print_latex_beam_particles(std::vector <HepMC::GenParticle*> beam_p){
        // there by default only two beam particles
        //
        // start vertices 
        std::cout   <<  "\\fmfleft{iv1,iv2}" << std::endl;
        // draw the lines: begin vertices are given by code above, the inner vertices are given by
        // "v"+endvertex(barcode)
        for (unsigned int i = 0; i< 2; i++){
            std::cout << "\\fmf{" 
                      << pdg_to_line(beam_p[i]->pdg_id()) 
                      << ",lab=$" << pdg_to_label(beam_p[i]->pdg_id()) << "$"  
                      <<"}{"
                      << "iv" << i+1 << ","
                      << "v" << -(beam_p[i]->end_vertex()->barcode())
                      << "}" << std::endl;
        }
}

void print_latex_interm_particles(std::vector <HepMC::GenParticle*> interm_p,std::map<HepMC::GenVertex *, int> level_map ,HepMC::GenEvent * evt){
    // see comments at print_latex_in_particles()
    //only lines
    for (unsigned int i = 0; i< interm_p.size(); i++){
        std::cout << "\\fmf{" 
                  << pdg_to_line(interm_p[i]->pdg_id()) ;
        if(level_map.count(interm_p[i]->production_vertex()) ){
            int level=level_map[interm_p[i]->production_vertex() ];
            double tension = pow(0.5,level*2);
        std::cout << ",tension=" << tension;
        }
        std::cout << ",lab=$" << pdg_to_label(interm_p[i]->pdg_id()) << "$"  
                  <<"}{";
        if (interm_p[i]->pdg_id()>0){
        std::cout << "v" << -(interm_p[i]->production_vertex()->barcode()) << ","
                  << "v" << -(interm_p[i]->end_vertex()->barcode()) ;
        }
        else{
        std::cout << "v" << -(interm_p[i]->end_vertex()->barcode()) << ","
                  << "v" << -(interm_p[i]->production_vertex()->barcode()) ;
        }
        std::cout << "}" << std::endl;
    }
    // let the beam particles end at fixed point
    if(evt->valid_beam_particles()){
        std::cout   << "\\fmfforce{(.2w,.9h)}{"
                    << "v" << -(evt->beam_particles().second->end_vertex()->barcode()) 
                    << "}" << std::endl;
        std::cout   << "\\fmfforce{(.2w,.1h)}{"
                    << "v" << -(evt->beam_particles().first->end_vertex()->barcode()) 
                    << "}" << std::endl;
    }        
    // let the interaction vertex be at a fixed point
    HepMC::GenVertex * iv = evt->signal_process_vertex();
    if (!iv) iv=find_interaction_vertex(evt);
    std::cout   << "\\fmfforce{(.3w,.5h)}{"
                << "v" << -(iv->barcode()) 
                << "}" << std::endl;
}

void print_latex_out_particles(std::vector <HepMC::GenParticle*> out_p, std::map<HepMC::GenVertex *, int> level_map, HepMC::GenEvent * evt){
        // see comments at print_latex_in_particles()
        //vertices
        std::cout << "\\fmfright{" ;
        for (unsigned int i=0; i< out_p.size(); i++){
            std::cout << "ov" << out_p[i]->barcode();
            if(i < (out_p.size()-1 )) std::cout << ","; 
        }
        std::cout << "}" << std::endl;
        //lines takes care of particles and antiparticles
        for (unsigned int i = 0; i< out_p.size(); i++){
            int level=level_map[out_p[i]->production_vertex() ];
            double tension = pow(0.5,level*2);
            std::cout << "\\fmf{" 
                      << pdg_to_line(out_p[i]->pdg_id()) 
                      << ",tension=" << tension
                      << ",lab=$" << pdg_to_label(out_p[i]->pdg_id()) << "$"  
                      <<"}{";
            if (out_p[i]->pdg_id()>0){
            std::cout << "v" << -(out_p[i]->production_vertex()->barcode()) << ","
                      << "ov" << out_p[i]->barcode() ;
            }
            else{
            std::cout << "ov" << out_p[i]->barcode() << ","
                      << "v" << -(out_p[i]->production_vertex()->barcode()) ;
            }
            std::cout << "}" << std::endl;
        }
}

void print_latex(HepMC::GenEvent* evt){
        std::vector <HepMC::GenParticle*> beam_p = get_beam_particles(evt);
        std::vector <HepMC::GenParticle*> interm_p = get_interm_particles(evt);
        std::pair< std::vector<HepMC::GenParticle*> , std::map<HepMC::GenVertex*, int> > out_p = get_sorted_out_particles(evt);
        std::cout << "\\fmfstraight" << std::endl ;
        print_latex_beam_particles(beam_p);
        print_latex_interm_particles(interm_p,out_p.second,evt);
        print_latex_out_particles(out_p.first,out_p.second,evt);
}


