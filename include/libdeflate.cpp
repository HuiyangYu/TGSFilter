#include <iostream>
#include <fstream>
#include <libdeflate.h>
#include <cstring>
#include <string>

#define OUTPUT_BUFFER_SIZE  5242880

using namespace std;

class DeflateCompress{
	public:
		DeflateCompress( ) {
			m_compressor = libdeflate_alloc_compressor(6);
		}

		~DeflateCompress() {
			libdeflate_free_compressor(m_compressor);
		}

		bool compressData(const void* input, size_t inputSize, uint8_t  ** m_outputBuffer, size_t & compressedSize, size_t & OUT_BUFFER_SIZE,int &  Thread) 
		{
			size_t bound = libdeflate_gzip_compress_bound(m_compressor, inputSize);

			if (OUT_BUFFER_SIZE<bound) {
				OUT_BUFFER_SIZE=bound;
				delete [] m_outputBuffer[Thread];
				m_outputBuffer[Thread]=new uint8_t[OUT_BUFFER_SIZE];
			}

			compressedSize = libdeflate_gzip_compress(m_compressor, input, inputSize, m_outputBuffer[Thread], OUT_BUFFER_SIZE);
			if (compressedSize == 0) {
				return false;
			}
			else {
				return true;
			}
		}

	private:
		libdeflate_compressor* m_compressor;
};

