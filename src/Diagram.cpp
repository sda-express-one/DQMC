#include "../include/Diagram.hpp"

//thread_local std::mt19937 Diagram::gen;

thread_local pcg32 Diagram::gen;

// constructor definition
Diagram::Diagram(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_int_max, int ph_ext_max, int data_type) : _N_diags(N_diags), _tau_max(tau_max),
    _chem_potential(chem_potential), _order_int_max(returnEven(order_int_max)), _ph_ext_max(ph_ext_max), _data_type(data_type) {

    // assign momentum values
    _kx = kx;
    _ky = ky;
    _kz = kz;

    if(_data_type == _data_type_array[0]){
        _vertices = new Vertex[_order_int_max + 2*_ph_ext_max + 2];
        _propagators = new Propagator[_order_int_max + 2*_ph_ext_max + 1];
    }
    else if(_data_type == _data_type_array[1]){
        _nodes = new FullVertexNode[_order_int_max + 2*_ph_ext_max + 2];

        // initialize array of all possible phonon vertices
        for(int i=0; i<_order_int_max + 2*_ph_ext_max + 2; i++){
            // phonon vertex specs
            _nodes[i].tau = 0;
            _nodes[i].tau_next = 0;
            _nodes[i].type = 0;
            _nodes[i].index = -1;

            // outgoing electron propagator specs
            _nodes[i].k[0] = _kx;
            _nodes[i].k[1] = _ky;
            _nodes[i].k[2] = _kz;

            // phonon momentum specs
            _nodes[i].w[0] = 0;
            _nodes[i].w[1] = 0;
            _nodes[i].w[2] = 0;

            if(i != 0){_nodes[i].prev = &_nodes[i-1];}
            if(i != _order_int_max + 2*_ph_ext_max + 1){_nodes[i].next = &_nodes[i+1];}
        }
        _free_list = &_nodes[0];
        _nodes[0].prev = nullptr;
        _nodes[_order_int_max + 2*_ph_ext_max + 1].next = nullptr;
    
        // collect and allocate first node of the diagram list
        _pointer_one = _free_list;
        _free_list = _free_list->next;
        _pointer_one->prev = nullptr;

        _pointer_one->next = nullptr;
        _diagram_list = _pointer_one;
        _pointer_one = nullptr;

        // fix starting node
        _head = _diagram_list;
        _tail = _diagram_list;


        //collect and allocate second node of the diagram list
        insertNode(_head);
        _tail = _diagram_list->next;
        _tail->tau = _tau_max/100;
        _tail->prev->tau_next = _tau_max/100;

        _internal_used = new FullVertexNodeIndicator[_order_int_max];
        _external_used = new FullVertexNodeIndicator[2*_ph_ext_max];

        for(int i = 0; i < std::max(_order_int_max, 2*_ph_ext_max); ++i){
            if(i < _order_int_max){_internal_used[i].position = i;}
            if(i < 2*_ph_ext_max){_external_used[i].position = i;}
        }
    }
};

Diagram::Diagram(FullVertexNode * nodes, int current_order, unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_int_max, int ph_ext_max, int data_type) : _N_diags(N_diags), _tau_max(tau_max),
    _chem_potential(chem_potential), _order_int_max(returnEven(order_int_max)), _ph_ext_max(ph_ext_max), _data_type(data_type) {

    // assign momentum values
    _kx = kx;
    _ky = ky;
    _kz = kz;

    _nodes = new FullVertexNode[_order_int_max + 2*_ph_ext_max + 2];
    _diagram_list = &_nodes[0];
    _head = &_nodes[0];
    _tail = &_nodes[1];
    for(int i = 0; i < _order_int_max + 2*_ph_ext_max + 2; ++i){
        if(i < current_order + 2){
            _nodes[i] = nodes[i];

            if(i == 0){
                _nodes[i].prev = nullptr;
                _nodes[i].next = &_nodes[2];
            }
            else if(i == 1){
                _nodes[i].prev = &_nodes[current_order+1];
                _nodes[i].next = nullptr;
            }
            else if(i == 2){
                _nodes[i].prev = &_nodes[0];
                _nodes[i].next = &_nodes[i+1];
            }
            else if(i == current_order + 1){
                _nodes[i].prev = &_nodes[i-1];
                _nodes[i].next = &_nodes[1];
            }
            else{
                _nodes[i].prev = &_nodes[i-1];
                _nodes[i].next = &_nodes[i+1];
            }
        }
        else{
            if(i == current_order + 2){
                _free_list = &_nodes[i];
                _nodes[i].prev = nullptr;
                _nodes[i].next = &_nodes[i+1];
            }
            else if(i == _order_int_max + 2*_ph_ext_max +1){
                _nodes[i].prev = &_nodes[i-1];
                _nodes[i].next = nullptr;
            }
            else{
                _nodes[i].prev = &_nodes[i-1];
                _nodes[i].next = &_nodes[i+1];
            }
        }
    }

    _internal_used = new FullVertexNodeIndicator[_order_int_max];
    _external_used = new FullVertexNodeIndicator[2*_ph_ext_max];


    // recollect links between arrays and node list (initialize support arrays)
    int j = 0;
    int internal_index = 0;
    int external_index = 0;
    bool found = false;
    // ugly double loop, if there is any other possibility please fix
    for(int i = 0; i < _order_int_max + 2*_ph_ext_max + 2; ++i){

        if(i < _order_int_max){_internal_used[i].position = i;}
        if(i < 2*_ph_ext_max){_external_used[i].position = i;}

        if(_nodes[i].type == +1){
            found = false;
            _internal_used[internal_index].linked = &_nodes[i];
            _internal_used[internal_index].used = true;
            _internal_used[internal_index].position = internal_index;
            j = 0;
            while(!found){
                if(_nodes[j].type == -1 && isEqual(_nodes[i].w[0], _nodes[j].w[0])){
                    found = true;
                    _internal_used[internal_index].conjugated = &_internal_used[internal_index + 1];
                    ++internal_index;
                    _internal_used[internal_index].linked = &_nodes[j];
                    _internal_used[internal_index].used = true;
                    _internal_used[internal_index].position = internal_index;
                    ++internal_index;
                }
                ++j;
            }
        }
        else if(_nodes[i].type == -2){
            found = false;
            _external_used[external_index].linked = &_nodes[i];
            _external_used[external_index].used = true;
            _external_used[external_index].position = external_index;
            j = 0;
            while(!found){
                if(_nodes[j].type == +2 && isEqual(_nodes[i].w[0], _nodes[j].w[0])){
                    found = true;
                    _external_used[external_index].conjugated = &_external_used[external_index + 1];
                    ++external_index;
                    _external_used[external_index].linked = &_nodes[j];
                    _external_used[external_index].used = true;
                    _external_used[external_index].position = external_index;
                    ++external_index;
                }
                ++j;
            }
        }
    }
    nodes = nullptr;
};

void Diagram::getNodes(FullVertexNode * nodes, int size){
    // check if array given in input is long enough
    if(size < _order_int_max + 2*_ph_ext_max + 2){return;}

    _helper = _head;
    int i = 2;
    while(_helper != nullptr){
        if(_helper == _head){
            nodes[0].electronic_band = _helper->electronic_band;
            nodes[0].index = _helper->index;
            nodes[0].tau = _helper->tau;
            nodes[0].tau_next = _helper->tau_next;
            nodes[0].type = _helper->type;
            nodes[0].k[0] = _helper->k[0];
            nodes[0].k[1] = _helper->k[1];
            nodes[0].k[2] = _helper->k[2];
            nodes[0].w[0] = _helper->w[0];
            nodes[0].w[1] = _helper->w[1];
            nodes[0].w[2] = _helper->w[2];
            nodes[0].prev = nullptr;
            nodes[0].next = &nodes[2];
        }
        else if(_helper == _tail){
            nodes[1].electronic_band = _helper->electronic_band;
            nodes[1].index = _helper->index;
            nodes[1].tau = _helper->tau;
            nodes[1].tau_next = 0;
            nodes[1].type = _helper->type;
            nodes[1].k[0] = _helper->k[0];
            nodes[1].k[1] = _helper->k[1];
            nodes[1].k[2] = _helper->k[2];
            nodes[1].w[0] = _helper->w[0];
            nodes[1].w[1] = _helper->w[1];
            nodes[1].w[2] = _helper->w[2];
            nodes[1].prev = &nodes[i-1];
            nodes[1].next = nullptr;
        }
        else{
            if(_helper->prev == _head){
                nodes[i].electronic_band = _helper->electronic_band;
                nodes[i].index = _helper->index;
                nodes[i].tau = _helper->tau;
                nodes[i].tau_next = _helper->tau_next;
                nodes[i].type = _helper->type;
                nodes[i].k[0] = _helper->k[0];
                nodes[i].k[1] = _helper->k[1];
                nodes[i].k[2] = _helper->k[2];
                nodes[i].w[0] = _helper->w[0];
                nodes[i].w[1] = _helper->w[1];
                nodes[i].w[2] = _helper->w[2];
                nodes[i].prev = &nodes[0];
                nodes[i].next = &nodes[i+1];
            }
            else if(_helper->next == _tail){
                nodes[i].electronic_band = _helper->electronic_band;
                nodes[i].index = _helper->index;
                nodes[i].tau = _helper->tau;
                nodes[i].tau_next = _helper->tau_next;
                nodes[i].type = _helper->type;
                nodes[i].k[0] = _helper->k[0];
                nodes[i].k[1] = _helper->k[1];
                nodes[i].k[2] = _helper->k[2];
                nodes[i].w[0] = _helper->w[0];
                nodes[i].w[1] = _helper->w[1];
                nodes[i].w[2] = _helper->w[2];
                nodes[i].prev = &nodes[i-1];
                nodes[i].next = &nodes[1];
            }
            else{
                nodes[i].electronic_band = _helper->electronic_band;
                nodes[i].index = _helper->index;
                nodes[i].tau = _helper->tau;
                nodes[i].tau_next = _helper->tau_next;
                nodes[i].type = _helper->type;
                nodes[i].k[0] = _helper->k[0];
                nodes[i].k[1] = _helper->k[1];
                nodes[i].k[2] = _helper->k[2];
                nodes[i].w[0] = _helper->w[0];
                nodes[i].w[1] = _helper->w[1];
                nodes[i].w[2] = _helper->w[2];
                nodes[i].prev = &nodes[i-1];
                nodes[i].next = &nodes[i+1];
            }
            ++i;
        }
        _helper = _helper->next;
    }
    nodes = nullptr;
};

void Diagram::insertNode(FullVertexNode * node_pointer){
    if(node_pointer == nullptr){return;}
    _helper = _free_list;
    _free_list = _free_list->next;
    if(_free_list != nullptr){_free_list->prev = nullptr;}
    _helper->prev = node_pointer;
    _helper->next = node_pointer->next;
    if(_helper->prev != nullptr){_helper->prev->next = _helper;}
    if(_helper->next != nullptr){_helper->next->prev = _helper;}
    _helper = nullptr;
};

void Diagram::deleteNode(FullVertexNode *& node_pointer){
    if(node_pointer == nullptr){return;}
    _helper = node_pointer;
    node_pointer = node_pointer->prev;
    node_pointer->next = _helper->next;
    _helper->next->prev = node_pointer;
    _helper->prev = nullptr;
    _helper->next = _free_list;
    if(_free_list != nullptr){
        _free_list->prev = _helper;
        _free_list = _free_list->prev;
    }
    else{
        _free_list = _helper;
    }
    _helper = nullptr;
};

// setters
void Diagram::setRelaxSteps(unsigned long long int relax_steps){_N_relax_steps = relax_steps;};

void Diagram::setAutcorrSteps(unsigned long long int autocorr_steps){_N_autocorr_steps = autocorr_steps;};