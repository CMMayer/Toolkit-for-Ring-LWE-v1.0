#include <iostream>

#include "Main_Classes/RingLweCryptographyField.h"
#include "Examples/CompactPublicKeyCryptosystem.h"
#include "Examples/DualStyleCryptosystem.h"
#include "Math util.h"

typedef RLWE_Toolkit::Main_Classes::RingLweCryptographyElement::Basis Basis;


/* 
A small main function that provides 4 examples. Example 1 and Example 2 use the implemented cryptosystems
to encrypt and decrypt the text "Bob,Alice". Example 3 and Example 4 provide small arithmetic examples
that can also be done by hand. This demonstrates that the computations are correct. 
*/
int main() {
	using namespace std;
	using namespace RLWE_Toolkit::Math_util;
	while (true)
	{
		cout << "---------------------------------------------------------------------------------------------- \n";
		cout << "Choose example by entering a number: \n"
			<< "0 = exit \n"
			<< "1 = compact public key cryptosystem \n"
			<< "2 = dual style cryptosystem\n"
			<< "3 = arithmetic 1\n"
			<< "4 = arithmetic 2\n"
			<< endl;

		int which_main = -1;
		cin >> which_main;
		cout << endl;

		pos_int m = 135;
		pos_int n = eulerTotient(m);
		pos_int p = 2;
		pos_int q = getPrimeMod1(m, 2000);

		switch (which_main)
		{
		case 1:
		{
			cout << "Example: compact public key cryptosystem \n\n";
			cout << "initializing...\n\n";

			CompactPublicKeyCryptosystem CPKC = CompactPublicKeyCryptosystem(m, p, q);

			char text[10] = { 'B', 'o', 'b', ',', 'A', 'l', 'i', 'c', 'e', '\0' };
			cout << "Plaintext: " << text << endl << endl;

			RingElement message = RingElement(CPKC.field_, Basis::BASIS_POWERFUL, toBitVec(text, n));
			cout << "Message in R_2: \n\n" << message;

			elmt_pair ciphertext = CPKC.enc(message);
			cout << "Message successfully encrypted! \n";

			char* enc_text = toCharArray(ciphertext.second->getCoordinates() % p);
			cout << "Last sequence of encrypted text: " << enc_text << endl;

			cout << "Decrypting message... \n \n";
			RingElement dec_message = CPKC.dec(ciphertext);
			char* dec_text = toCharArray(dec_message.getCoordinates());
			cout << "Decrypted message: " << dec_text << endl << endl;
		}
		break;

		case 2:
		{
			cout << "Example: dual style cryptosystem \n\n";
			cout << "initializing...\n\n";

			int l = 2;
			DualStyleCryptosystem DSC = DualStyleCryptosystem(m, p, q, l);

			char text[10] = { 'B', 'o', 'b', ',', 'A', 'l', 'i', 'c', 'e', '\0' };
			cout << "Plaintext: " << text << endl << endl;

			RingElement message = RingElement(DSC.field_, Basis::BASIS_POWERFUL, toBitVec(text, n));
			cout << "Message in R_2: \n\n" << message;

			elmt_vec ciphertext = DSC.enc(message);
			cout << "Message successfully encrypted! \n";

			char* enc_text = toCharArray(ciphertext[l - 1].getCoordinates() % p);
			cout << "Last sequence of encrypted text: " << enc_text << endl;

			cout << "Decrypting message... \n \n";
			RingElement dec_message = DSC.dec(ciphertext);
			char* dec_text = toCharArray(dec_message.getCoordinates());
			cout << "Decrypted message: " << dec_text << endl << endl;

		}
		break;

		case 3:
		{
			cout << "Example: arithmetic 1 \n\n";

			cout << "Initializing cyclotomic field K with m = 9. Therefore p = 3 is the only prime \n"
				<< "dividing m and the powerful basis of K is (XI_9^0, ..., XI_9^5). This small basis allows us, \n"
				<< "to create examples that can also be computed by hand. \n\n";
			m = 9;
			q = getPrimeMod1(m, m);

			shared_ptr<RingLweCryptographyField> field = make_shared<RingLweCryptographyField>(m, q);

			cout << "Example 1: \n"
				<< "Multiply the elements given by (1,0,-1,0,0,0) and (0,1,1,1,0,0). The expected coordinate vector is (0,1,1,0,-1,-1)\n\n";
			coordinate_vec c1 = { 1, 0, -1, 0, 0, 0 };
			coordinate_vec c2 = { 0, 1, 1, 1, 0, 0 };
			RingElement test1 = RingElement(field, Basis::BASIS_POWERFUL, c1);
			RingElement test2 = RingElement(field, Basis::BASIS_POWERFUL, c2);
			cout << "The computed element is:\n\n" << test1*test2;

			cout << "Example 2: \n"
				<< "Square the element given by (0,1,1,0,0,0). The expected coordinate vector is (0,0,1,2,1,0)\n\n";
			coordinate_vec c3 = { 0, 1, 1, 0, 0, 0 };
			RingElement test3 = RingElement(field, Basis::BASIS_POWERFUL, c3);
			cout << "The computed element is:\n\n" << test3*test3;

			cout << "Example 3: \n"
				<< "Multiply the elements given by (0,2,-5,1,0,0) and (3,4,-7,0,0,0). The expected coordinate vector is (0,6,-7,-31,39,-7)\n\n";
			coordinate_vec c4 = { 0, 2, -5, 1, 0, 0 };
			coordinate_vec c5 = { 3, 4, -7, 0, 0, 0 };
			RingElement test4 = RingElement(field, Basis::BASIS_POWERFUL, c4);
			RingElement test5 = RingElement(field, Basis::BASIS_POWERFUL, c5);
			cout << "The computed element is:\n\n" << test4*test5;
		}
		break;

		case 4:
		{
			cout << "Example: arithmetic 2 \n\n";

			cout << "Initializing cyclotomic field K with m = 9. Therefore p = 3 is the only prime \n"
				<< "dividing m and the powerful basis of K is (XI_9^0, ..., XI_9^5). An easy computation shows that \n"
				<< "the decoding basis for R dual is 1/9(1 - XI_9^6, XI_9^1 - XI_9^7, XI_9^2 - XI_9^8, XI_9^3 - XI_9^6, XI_9^4 - XI_9^7, XI_9^5 - XI_9^8) \n"
				<< "This small basis allows us to create examples that can also be computed by hand. \n\n";
			m = 9;
			q = getPrimeMod1(m, m);

			shared_ptr<RingLweCryptographyField> field = make_shared<RingLweCryptographyField>(m, q);

			cout << "Example 1: \n"
				<< "Add the elements given by (0,0,1,0,0,0) in R dual and (0,1,0,0,0,0) in R. The element is thus given by XI_9 + 1/9(XI_9^2 - XI_9^8)\n\n";
			coordinate_vec c1 = { 0, 0, 1, 0, 0, 0 };
			coordinate_vec c2 = { 0, 1, 0, 0, 0, 0 };
			RingElement test1 = RingElement(field, Basis::BASIS_DECODING, c1);
			RingElement test2 = RingElement(field, Basis::BASIS_POWERFUL, c2);
			cout << "The computed element is:\n\n" << test1 + test2
				<< "A quick computation shows that (0,6,1,0,-3,0) is indeed the coordinate vector in the decoding basis.\n\n";

			cout << "Example 2: \n"
				<< "Taking the same elements and multiply them should give the coordinate vector (-1,0,0,1,0,0) \n";
			cout << "The computed element is:\n\n" << test1 * test2;
		}
			break;

		case 0:
			return 0;

		default:
			return 0;

		break;
		}

	}

}