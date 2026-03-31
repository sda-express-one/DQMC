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

Diagram::Diagram(FullVertexNode * nodes, FullVertexNodeIndicator * internal_used, FullVertexNodeIndicator * external_used, 
    int current_order, unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_int_max, int ph_ext_max, int data_type) : _N_diags(N_diags), _tau_max(tau_max),
    _chem_potential(chem_potential), _order_int_max(returnEven(order_int_max)), _ph_ext_max(ph_ext_max), _data_type(data_type) {

    // assign momentum values
    _kx = kx;
    _ky = ky;
    _kz = kz;

    _internal_used = new FullVertexNodeIndicator[_order_int_max];
    _external_used = new FullVertexNodeIndicator[2*_ph_ext_max];

    _nodes = new FullVertexNode[_order_int_max + 2*_ph_ext_max + 2];
    _diagram_list = &_nodes[0];
    _head = &_nodes[0];
    _tail = &_nodes[1];
    
    int j = 0, k = 0;
    int counter_int = 0, counter_ext = 0;
    FullVertexNode * second_helper = nullptr;
    FullVertexNode * third_helper = nullptr;

    if(current_order == 0){
        _nodes[0] = nodes[0];
        _nodes[1] = nodes[1];
        _nodes[0].prev = nullptr;
        _nodes[0].next = &_nodes[1];
        _nodes[1].prev = &_nodes[0];
        _nodes[1].next = nullptr;

        _internal_used[0].position = 0;
        _external_used[0].position = 0;
        _internal_used[1].position = 1;
        _external_used[1].position = 1;

        for(int i = 2; i < _order_int_max + 2*_ph_ext_max + 2; ++i){
            if(i < _order_int_max){_internal_used[i].position = i;}
            if(i < 2*_ph_ext_max){_external_used[i].position = i;}

            if(i == 2){
                _free_list = &_nodes[i];
                _nodes[i].prev = nullptr;
                _nodes[i].next = &_nodes[i+1];
            }
            else if(i == _order_int_max + 2*_ph_ext_max + 1){
                _nodes[i].prev = &_nodes[i-1];
                _nodes[i].next = nullptr;
            }   
            else{
                _nodes[i].prev = &_nodes[i-1];
                _nodes[i].next = &_nodes[i+1];
            }
        }
    }
    else{
        for(int i = 0; i < _order_int_max + 2*_ph_ext_max + 2; ++i){
            if(i < _order_int_max){_internal_used[i].position = i;}
            if(i < 2*_ph_ext_max){_external_used[i].position = i;}

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

                
                if(_nodes[i].type == +1){
                    second_helper = &nodes[i];
                    j = 0;
                    k = 0;
                    
                    while(internal_used[j].linked != second_helper){++j;}
                    _internal_used[counter_int].linked = &_nodes[i];
                    _internal_used[counter_int].used = true;
                    ++counter_int;

                    third_helper = internal_used[j].conjugated->linked;
                    while(&nodes[k] != third_helper){++k;}
                    _internal_used[counter_int].linked = &_nodes[k];
                    _internal_used[counter_int].used = true;

                    _internal_used[counter_int].conjugated = &_internal_used[counter_int - 1];
                    _internal_used[counter_int - 1].conjugated = &_internal_used[counter_int];
                    ++counter_int;
                }
                else if(_nodes[i].type == -2){
                    second_helper = &nodes[i];
                    j = 0;
                    k = 0;

                    while(external_used[j].linked != second_helper){++j;}
                    _external_used[counter_ext].linked = &_nodes[i];
                    _external_used[counter_ext].used = true;
                    ++counter_ext;

                    third_helper = external_used[j].conjugated->linked;
                    while(&nodes[k] != third_helper){++k;}
                    _external_used[counter_ext].linked = &_nodes[k];
                    _external_used[counter_ext].used = true;

                    _external_used[counter_ext].conjugated = &_external_used[counter_ext - 1];
                    _external_used[counter_ext - 1].conjugated = &_external_used[counter_ext];
                    ++counter_ext;
                }

                second_helper = nullptr;
                third_helper = nullptr;
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
    }
    nodes = nullptr;
    internal_used = nullptr;
    external_used = nullptr;
};

void Diagram::getNodes(FullVertexNode * nodes, FullVertexNodeIndicator * internal_used, FullVertexNodeIndicator * external_used, int size){
    // check if array given in input is long enough
    if(size < _order_int_max + 2*_ph_ext_max + 2){return;}

    for(int i = 0; i < std::max(_order_int_max, 2*_ph_ext_max); ++i){
        if(i < _order_int_max){internal_used[i].position = i;}
        if(i < 2*_ph_ext_max){external_used[i].position = i;}
    }

    _helper = _head;
    FullVertexNode * second_helper = nullptr;
    FullVertexNode * third_helper = nullptr;
    int i = 2;
    int j = 0, k = 0;
    int counter_int = 0, counter_ext = 0;
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

            if(_helper->type == +1){
                j = 0;
                while(_internal_used[j].linked != _helper){++j;}
                second_helper = _internal_used[j].conjugated->linked;
                internal_used[counter_int].linked = &nodes[i];
                internal_used[counter_int].used = true;
                ++counter_int;

                k = 0;
                third_helper = _head;
                while(third_helper != second_helper && third_helper != nullptr){
                    third_helper = third_helper->next;
                    ++k;
                }
                internal_used[counter_int].linked = &nodes[k+1]; // to be checked
                internal_used[counter_int].used = true;
                internal_used[counter_int].conjugated = &internal_used[counter_int-1];
                internal_used[counter_int-1].conjugated = &internal_used[counter_int];
                ++counter_int;

                second_helper = nullptr;
                third_helper = nullptr;
            }
            else if(_helper->type == -2){
                j = 0;
                while(_external_used[j].linked != _helper){++j;}
                second_helper = _external_used[j].conjugated->linked;
                external_used[counter_ext].linked = &nodes[i];
                external_used[counter_ext].used = true;
                ++counter_ext;

                k = 0;
                third_helper = _head;
                while(third_helper != second_helper && third_helper != nullptr){
                    ++k;
                    third_helper = third_helper->next;
                }
                external_used[counter_ext].linked = &nodes[k+1];
                external_used[counter_ext].used = true;
                external_used[counter_ext].conjugated = &external_used[counter_ext-1];
                external_used[counter_ext-1].conjugated = &external_used[counter_ext];
                ++counter_ext;

                second_helper = nullptr;
                third_helper = nullptr;
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