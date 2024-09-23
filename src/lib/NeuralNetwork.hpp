#ifndef NN_H
#define NN_H

#include <vector>
#include <functional>
#include <string.h>
#include <fstream>
#include <iostream>
#include "Validate.hpp"

// Special thanks to the CS department - Rooms 13 & 14
const unsigned file_signature = 0x43531314;

const std::function<double(double)> af_linear = [](double value_) {return value_;};

class Layer {
    
    // number of nodes in this layer
    const unsigned node_count;
    // number of input nodes
    const unsigned input_node_count;

    public:
    std::vector<double> bias;
    std::vector<double> weight;

    std::function<double(double)> activation_function;

    Layer(unsigned node_count_, unsigned input_count_, std::function<double(double)> af_ = af_linear) :
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
            output[i] = activation_function(output[i] += bias[i]);
        }
        return output;
    }
};

struct NeuralNetwork {

    std::vector<Layer> layers;

    NeuralNetwork(std::vector<unsigned> shape_) {
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
        for (auto layer : layers) {
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