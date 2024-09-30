#ifndef NN_H
#define NN_H

#include <vector>
#include <cmath>
#include <functional>
#include <string.h>
#include <fstream>
#include <iostream>
#include "Validate.hpp"

// Special thanks to the CS department - Rooms 13 & 14
const unsigned file_signature = 0x43531314;

const double e_ = 2.718281828459045;

// NOTE - These functions do not need to be pointers, but it seems to cause undefined behavior when they're not pointers

// linear activation function
// range: x ∈ R
const auto af_linear = new std::function([](double value_) {return value_;});
// sigmoid activation function
// range: -1 < x < 1
const auto af_sig = new std::function([](double value_) {
    // https://en.wikipedia.org/wiki/Sigmoid_function
    double e_pow_neg_x = pow(e_,-value_);
    return (1.0-e_pow_neg_x)/(1.0+e_pow_neg_x);
});
// sigmoid activation function
// range: 0 < x < 1
const auto af_sig01 = new std::function([](double value_) {
    // https://en.wikipedia.org/wiki/Sigmoid_function
    return 1.0/(1.0+pow(e_,-value_));
});

class Layer {
    
    // number of nodes in this layer
    const unsigned node_count;
    // number of input nodes
    const unsigned input_node_count;

    public:
    std::vector<double> bias;
    std::vector<double> weight;

    // pointer to the activation function
    std::function<double(double)>* activation_function;

    Layer(unsigned node_count_, unsigned input_count_, std::function<double(double)>* af_ = af_linear) :
    node_count(node_count_), input_node_count(input_count_), bias(node_count_, 0.0), weight(node_count_ * input_count_, 0.0), activation_function(af_)
    {}
    /**
     * \brief Calculate the output of the layer
     * \param input_ Incoming input
     * \return The layer output
    */
    std::vector<double> calculate(std::vector<double> input_) {
        std::vector<double> output(node_count, 0.0);
        unsigned inp_i = 0; // input index
        unsigned out_i = 0; // output index
        for (auto w : weight) {
            output[out_i] += w*input_[inp_i];
            out_i++;
            // if the output index has reached the end
            // start back at the first node and increment input
            if (out_i == node_count) {
                out_i = 0;
                inp_i++;
            }
        }
        for (unsigned i = 0; i < node_count; i++) {
            // apply bias and activation function
            output[i] = (*activation_function)(output[i] += bias[i]);
        }
        return output;
    }
};

struct NeuralNetwork {

    std::vector<Layer> layers;
    unsigned input_count;
    unsigned output_count;

    NeuralNetwork(std::vector<unsigned> shape_) : input_count(shape_[0]), output_count(shape_[shape_.size()-1]) {
        unsigned incoming_input_count = shape_[0];
        for (unsigned i = 1; i < shape_.size(); i++) {
            Layer layer(shape_[i], incoming_input_count);
            incoming_input_count = shape_[i];
            layers.push_back(layer);
        }
    }
    /**
     * \brief Calculate the output of the network
     * \param input_ Incoming input
     * \return The network output
    */
    std::vector<double>calculate(std::vector<double>& input_) {
        std::vector<double> output = input_;
        for (auto& layer : layers) {
            output = layer.calculate(output);
        }
        return output;
    }
    /**
     * \brief Save network values to a file
     * \param path_ File path
    */
    void save_to_file(char* path_) {
        std::ofstream file(path_, std::ios::binary);
        file.write(reinterpret_cast<const char*>(&file_signature), sizeof(unsigned));
        // https://www.reddit.com/r/cpp_questions/comments/7i50i3/what_is_an_elegant_modern_way_to_quickly_read_and/
        for (auto& layer : layers) {
            file.write(reinterpret_cast<const char*>(layer.weight.data()), layer.weight.size() * sizeof(double));
            file.write(reinterpret_cast<const char*>(layer.bias.data()), layer.bias.size() * sizeof(double));
        }
        file.close();
    }
    /**
     * \brief Load network values from a file
     * \note Loading values from a network of a different size is undefined behavior
     * \param path_ File path
    */
    void load_from_file(char* path_) {
        std::ifstream file(path_, std::ios::binary);
        unsigned sig;
        file.read(reinterpret_cast<char*>(&sig), sizeof(unsigned));
        if (sig != file_signature) {
            send_error_message(path_, error_type::FORMAT_ERROR, ERROR);
            throw;
        }
        // https://www.reddit.com/r/cpp_questions/comments/7i50i3/what_is_an_elegant_modern_way_to_quickly_read_and/
        for (auto& layer : layers) {
            file.read(reinterpret_cast<char*>(layer.weight.data()), layer.weight.size() * sizeof(double));
            file.read(reinterpret_cast<char*>(layer.bias.data()), layer.bias.size() * sizeof(double));
        }
        file.close();
    }
};

#endif