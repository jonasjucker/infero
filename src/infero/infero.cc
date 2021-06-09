/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <memory>

#include "eckit/log/Log.h"
#include "eckit/option/CmdArgs.h"
#include "eckit/option/SimpleOption.h"
#include "eckit/runtime/Main.h"
#include "eckit/serialisation/FileStream.h"
#include "eckit/config/LocalConfiguration.h"

#include "infero/inference_models/InferenceModel.h"
#include "infero/infero_utils.h"


using namespace eckit;
using namespace eckit::option;
using namespace eckit::linalg;

using namespace infero;
using namespace infero::utils;


void usage(const std::string&) {

    Log::info() << std::endl
                << "-------------------------------" << std::endl
                << "Machine Learning execution tool" << std::endl
                << "-------------------------------" << std::endl
                << std::endl
                << "Reads a ML model + input data and "
                   "generates a prediction."
                << std::endl;
}


int main(int argc, char** argv) {

    Main::initialise(argc, argv);
    std::vector<Option*> options;

    options.push_back(new SimpleOption<std::string>("input", "Path to input file"));
    options.push_back(new SimpleOption<std::string>("output", "Path to output file"));
    options.push_back(new SimpleOption<std::string>("model", "Path to ML model"));
    options.push_back(new SimpleOption<std::string>("engine", "ML engine [onnx, tflite, trt]"));
    options.push_back(new SimpleOption<std::string>("ref_path", "Path to Reference prediction"));
    options.push_back(new SimpleOption<double>("threshold", "Verification threshold"));

    CmdArgs args(&usage, options, 0, 0, true);

    // Input path
    std::string input_path  = args.getString("input", "data.npy");
    std::string model_path  = args.getString("model", "model.onnx");
    std::string engine_type = args.getString("engine", "onnx");
    std::string output_path = args.getString("output", "out.csv");
    std::string ref_path    = args.getString("ref_path", "");
    float threshold         = args.getFloat("threshold", 0.001f);

    // Model configuration from CL
    LocalConfiguration local;
    local.set("path", model_path);

    // Input data
    TensorFloat* inputT = tensor_from_file<float>(input_path);
    TensorFloat predT;

    // Inference model
    auto engine = InferenceModel::open(engine_type, local);
    std::cout << *engine << std::endl;

    // Run inference
    engine->infer(*inputT, predT);

    // Save output tensor to file
    if (args.has("output")) {
        tensor_to_file<float>(predT, output_path);
    }

    // Compare against ref values
    if (args.has("ref_path")) {

        TensorFloat* refT = tensor_from_file<float>(ref_path);
        float err = compare_tensors<float>(predT, *refT, TensorErrorType::MSE);
        Log::info() << "MSE error: " << err << std::endl;
        Log::info() << "threshold: " << threshold << std::endl;

        delete refT;
        delete inputT;
        return !(err < threshold);
    }

    delete inputT;

    return EXIT_SUCCESS;
}
