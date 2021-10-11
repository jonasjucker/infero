/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
typedef void* infero_model_handle;
void infero_initialise(int argc, char** argv);
infero_model_handle infero_create_handle_from_yaml_str(char str[]);
infero_model_handle infero_create_handle_from_yaml_file(char path[]);
void infero_open_handle(infero_model_handle);
void infero_close_handle(infero_model_handle);
void infero_delete_handle(infero_model_handle);
void infero_inference_double(infero_model_handle h, double data1[], int rank1, int shape1[], double data2[], int rank2, int shape2[]);
void infero_inference_float(infero_model_handle h, float data1[], int rank1, int shape1[], float data2[], int rank2, int shape2[]);
void infero_finalise();