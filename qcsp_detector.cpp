/* MyMEXFunction
 * c = MyMEXFunction(a,b);
 * Adds offset argument a to each element of double array b and
 * returns the modified array c.
 */

#include "includes/CDetectionStateGeneric.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <complex>
#include <stdexcept>
#include <string>

using namespace matlab::data;
using matlab::mex::ArgumentList;

#include "includes/CDetectorSerialGeneric.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

struct scalarArguments {
    std::vector<float> pn;
    int                N;
    int                p_omega;
    float              threshold;
    int                step_denominator;
    int                step_numerator;
    bool               normed;
};

class StateVector {
private:
    const uint16_t        _p_omega;
    size_t                _size;
    std::vector<bool>     _frame_detected;
    std::vector<bool>     _max_found;
    std::vector<float>    _scores; // Column major linearized scores (x[i][j] <=> x[i*p_omega + j])
    std::vector<float>    _max_score;
    std::vector<uint64_t> _chip_since_last_det;
    std::vector<uint64_t> _chip_from_max;
    std::vector<float>    _frequency_offset; // The frequency offset corresponding to the maximum score
    std::vector<uint32_t> _frequency_index;  // The index of the frequency offset in the frequency array

public:
    using State = QCSP::StandaloneDetector::DetectionStateGeneric;
    const uint16_t &              p_omega() const { return _p_omega; }
    const size_t &                size() const { return _size; }
    const std::vector<bool> &     frame_detected() const { return _frame_detected; }
    const std::vector<bool> &     max_found() const { return _max_found; }
    const std::vector<float> &    scores() const { return _scores; }
    const std::vector<float> &    max_score() const { return _max_score; }
    const std::vector<uint64_t> & chip_since_last_det() const { return _chip_since_last_det; }
    const std::vector<uint64_t> & chip_from_max() const { return _chip_from_max; }
    const std::vector<float> &    frequency_offset() const { return _frequency_offset; } // The frequency offset corresponding to the maximum score
    const std::vector<uint32_t> & frequency_index() const { return _frequency_index; }   // The index of the frequency offset in the frequency array

    void push_back(const State & s) {
        this->_frame_detected.push_back(s.frame_detected);
        this->_max_found.push_back(s.max_found);
        this->_max_score.push_back(s.max_score);
        this->_chip_since_last_det.push_back(s.chip_since_last_det);
        this->_chip_from_max.push_back(s.chip_from_max);
        this->_frequency_offset.push_back(s.frequency_offset);
        this->_frequency_index.push_back(s.frequency_index);
        for (int i = 0; i < this->p_omega(); i++) {
            this->_scores.push_back(s.scores[i]);
        }
        _size++;
    }

    void set(size_t idx, const State & s) {
        if (idx >= _size) {
            throw std::runtime_error(
                "Index " + std::to_string(idx) + " overcome the capacity (" + std::to_string(this->size()) + ").");
        }
        this->_frame_detected[idx]      = s.frame_detected;
        this->_max_found[idx]           = s.max_found;
        this->_max_score[idx]           = s.max_score;
        this->_chip_since_last_det[idx] = s.chip_since_last_det;
        this->_chip_from_max[idx]       = s.chip_from_max;
        this->_frequency_offset[idx]    = s.frequency_offset;
        this->_frequency_index[idx]     = s.frequency_index;
        for (int i = 0; i < this->p_omega(); i++) {
            this->_scores[i + idx] = s.scores[i];
        }
    }

    StateVector() = delete;
    StateVector(uint16_t p_omega, size_t size = 0)
        : _p_omega(p_omega),
          _size(size),
          _frame_detected(size, false),
          _max_found(size, false),
          _scores(p_omega * size, 0.f),
          _max_score(size, 0.f),
          _chip_since_last_det(size, 0LU),
          _chip_from_max(size, 0),
          _frequency_offset(size, 0.f),
          _frequency_index(size, 0U) {}
    StateVector(const std::vector<State> & v)
        : _p_omega(v.at(0).p_omega),
          _size(v.size()),
          _frame_detected(this->_size, false),
          _max_found(this->_size, false),
          _scores(this->_p_omega * this->_size, 0.f),
          _max_score(this->_size, 0.f),
          _chip_since_last_det(this->_size, 0LU),
          _chip_from_max(this->_size, 0),
          _frequency_offset(this->_size, 0.f),
          _frequency_index(this->_size, 0U) {
        for (int i = 0; i < this->size(); i++) {
            this->set(i, v[i]);
        }
    }
};

class MexFunction : public matlab::mex::Function {
private:
    ArrayFactory f;

    StateVector * a(
        const std::vector<std::complex<float>> & inData,
        const scalarArguments &                  args) {
        QCSP::StandaloneDetector::CDetectorSerialGeneric detector(
            args.pn,
            args.N,
            args.p_omega,
            args.threshold,
            args.step_denominator,
            args.step_numerator,
            args.normed);

        StateVector * psv = new StateVector(args.p_omega);

        const float *      inRawPtr = reinterpret_cast<const float *>(inData.data());
        StateVector::State s(args.p_omega); // Initial state
        for (int i = 0; i < inData.size(); i++) {
            const int idx = i << 1;
            detector.process(inRawPtr[idx], inRawPtr[idx + 1], &s);
            psv->push_back(s);
        }

        return psv;
    }

public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        // checkArguments(outputs, inputs);

        const TypedArray<std::complex<float>>  cpxDataMat = std::move(inputs[0]);
        const std::vector<std::complex<float>> cpxData(cpxDataMat.begin(), cpxDataMat.end());

        const TypedArray<float>  pnMat = std::move(inputs[1]);
        const std::vector<float> pn(pnMat.begin(), pnMat.end());

        const size_t nInputs = inputs.size();

        scalarArguments args {pn, 60, 1, 140.f, 1, 1, true};

        if (nInputs > 2)
            args.N = int(inputs[2][0]);
        if (nInputs > 3)
            args.p_omega = int(inputs[3][0]);
        if (nInputs > 4)
            args.threshold = float(inputs[4][0]);
        if (nInputs > 5)
            args.step_denominator = int(inputs[5][0]);
        if (nInputs > 6)
            args.step_numerator = int(inputs[6][0]);
        if (nInputs > 7)
            args.normed = bool(inputs[7][0]);

        // TypedArray<double> doubleArray = std::move(inputs[1]);

        StateVector * psv = a(cpxData, args);

        const size_t nSize = psv->size();

        const size_t nOutputs = outputs.size();
        if (nOutputs > 0) {
            outputs[0] = f.createArray<float>(
                {1, nSize},
                psv->max_score().data(),
                psv->max_score().data() + nSize);
        }
        if (nOutputs > 1) {
            outputs[1] = f.createArray<float>(
                {1, nSize},
                psv->frequency_offset().data(),
                psv->frequency_offset().data() + nSize);
        }
        if (nOutputs > 2) {
            outputs[2] = f.createArray<float>({unsigned(args.p_omega), nSize},
                                              psv->scores().data(),
                                              psv->scores().data() + psv->scores().size());
        }
        if (nOutputs > 3) {
            bool * ptr = new bool[nSize];
            std::copy(psv->frame_detected().begin(), psv->frame_detected().end(), ptr);
            outputs[3] = f.createArray<bool>(
                {1, nSize},
                ptr,
                ptr + nSize);
        }
        if (nOutputs > 4)
            outputs[4] = f.createArray<uint32_t>(
                {1, nSize},
                psv->frequency_index().data(),
                psv->frequency_index().data() + nSize);
        if (nOutputs > 5) {
            bool * ptr = new bool[nSize];
            std::copy(psv->max_found().begin(), psv->max_found().end(), ptr);
            outputs[5] = f.createArray<bool>(
                {1, nSize},
                ptr,
                ptr + nSize);
        }
        if (nOutputs > 6)
            outputs[6] = f.createArray<uint64_t>(
                {1, nSize},
                psv->chip_from_max().data(),
                psv->chip_from_max().data() + nSize);
        if (nOutputs > 7)
            outputs[7] = f.createArray<uint64_t>(
                {1, nSize},
                psv->chip_since_last_det().data(),
                psv->chip_since_last_det().data() + nSize);
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        // Get pointer to engine
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

        // Get array factory
        ArrayFactory factory;

        // Check offset argument: First input must be scalar double
        if (inputs[0].getType() != ArrayType::DOUBLE || inputs[0].getNumberOfElements() != 1) {
            matlabPtr->feval(u"error",
                             0,
                             std::vector<Array>({factory.createScalar("First input must be scalar double")}));
        }

        // Check array argument: Second input must be double array
        if (inputs[1].getType() != ArrayType::DOUBLE) {
            matlabPtr->feval(u"error",
                             0,
                             std::vector<Array>({factory.createScalar("Input must be double array")}));
        }
        // Check number of outputs
        if (outputs.size() > 1) {
            matlabPtr->feval(u"error",
                             0,
                             std::vector<Array>({factory.createScalar("Only one output is returned")}));
        }
    }
};