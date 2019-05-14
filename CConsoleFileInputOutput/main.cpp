#include <stdio.h>
#include <iostream>
#include "math.h"

int main(int argc, char **argv)
{
/*
//C input/output
    //printf:
    int a,b,c;        // Объявление переменных a,b,c
    a=5;
    b=6;
    c=9;
    printf("a=%d, b=%d, c=%d\n",a,b,c);
    //Точность задаётся следующим образом: %.n<код формата>.
    //Где n - число цифр после запятой, а <код формата> - один из кодов приведённых выше.
    float x,y,z;
    x=10.5;
    y=130.67;
    z=54;
    printf("Coordinates: x:%.2f, y:%.2f, z:%.2f\n", x, y, z);
    //
    printf("%5d",20);//output:   20
    printf("%05d",20);//output:00020
    //%% - %; /" - кавычки
    a=11;        // 11 в десятичной равно b в шестнадцатеричной
    //%d- decimal %X - hexadecimal
    printf("a-dec=%d, a-hex=%X",a,a);
    //
    char ch1,ch2,ch3;
    ch1='A';
    ch2='B';
    ch3='C';
    printf("%c%c%c",ch1,ch2,ch3);
    //
    char *str="My line.";
    printf("This is %s",str);
    //scanf:
    //При вводе строки с помощью функции scanf() (спецификатор формата %s),
    //строка вводиться до первого пробела!! т.е. если вы вводите строку
    //"Привет мир!" с использованием функции scanf()
    char str1[80];        // массив на 80 символов
    scanf("%s",str1);
    //то после ввода результирующая строка,
    //которая будет храниться в массиве str будет состоять из одного слова "Привет".
    //ФУНКЦИЯ ВВОДИТ СТРОКУ ДО ПЕРВОГО ПРОБЕЛА!

    //Если вы хотите вводить строки с пробелами, то используйте функцию
    //char *gets( char *buf );
    //С помощью функции gets() вы сможете вводить полноценные строки.
    //Функция gets() читает символы с клавиатуры до появления символа
    //новой строки (\n). Сам символ новой строки появляется, когда вы
    //нажимаете клавишу enter. Функция возвращает указатель на buf.
    //buf - буфер (память) для вводимой строки.
    char buffer2[100];       // массив (буфер) для вводимой строки
    gets(buffer2);            // вводим строку и нажимаем enter
    printf("%s",buffer2);    // вывод введённой строки на экран
    //Для ввода данных с помощью функции scanf(), ей в качестве параметров
    //нужно передавать адреса переменных, а не сами переменные
    int x1;
    printf("Type x1:");
    scanf("%d",&x1);
    printf("x1=%d",x1);
    //
    int age;
    printf("\nHow old are you:");
    scanf("%d",&age);
    printf("You are %d years old.", age);
    //
    int x2, y2;
    printf("\nCalculator:");
    scanf("%d+%d", &x2, &y2);
    printf("\n%d+%d=%d", x2, y2, x2+y2);
    //
    char name[5];
    printf("\nType your login (not more than 5 symbols):");
    scanf("%5s", name);
    printf("\nYou typed %s", name);
    //как нужно использовать множество поиска
    char bal;
    printf("Your mark 2,3,4,5:");
    scanf("%[2345]", &bal);
    printf("\nMark %c", bal);

 */

 /*
//------1.enter single value in c-------------------------------------------------//
    int c;
    printf( "Enter a value :");
    c = getchar( );
    printf( "\nYou entered: \n");
    putchar( c );
    //exit(0);
//--1.enter single value in c--and--2.enter string in c---do not work together--//
//--only one of then works--//
//----------2.enter string in c-----------------------------------------//
    char str[100];
    //printf("H1");
    printf( "\nEnter a string:\n");
    gets( str );
    //printf("H2");
    printf( "\nYou entered: \n");
    puts( str );
    //printf("H3");
    //exit(0);

//-------------------------------------------------------------------------------//
//Между тем использовать функцию gets категорически не рекомендуется, ввиду того,
//что она не контролирует выход за границу строки, что может произвести к ошибкам.
//Вместо нее используется функция fgets с тремя параметрами:
//    char * fgets(char * buffer, int size, FILE * stream);
//-------------------------------------------------------------------------------//
    char str3[102] = "";
    printf("Enter a string: ");
    fgets(str3, 102, stdin);
    printf("You entered: ");
    puts(str3);
    exit(0);
//-------------------------------------------------------------------------------//

*/

//C file input/output:
    FILE *myfile;
//it is written in the file:
//apples 10 23.4
//bananas 5 25.0
//bread 1 10.3
    FILE *file;
    struct food
    {
       char * name;
       unsigned qty;
       float price;
    };
    struct food shop[3];

    std::cout<<"HERE WE ARE"<<std::endl;
//----------------------------------------------------------//
//print 2 numbers to file
    int n1 = 16;
    float n2 = 3.141592654f;
    FILE *fp;
    fp = fopen( "data.txt", "w" );
    fprintf( fp, "  %d   %f", n1, n2 );//"  %d   %f" - format line
    fclose( fp );
//read 2 numbers from file
    int rn1;
    float rn2;
    fp = fopen( "data.txt", "r" );
    fscanf( fp, "%d %f", &rn1, &rn2 );
    printf( "%d %f", rn1, rn2 );
    fclose( fp );
//----------------------------------------------------------//
//https://en.wikibooks.org/wiki/A_Little_C_Primer/C_File-IO_Through_Library_Functions
//-------write and read:-----------------------------------------------//
    int ctr, ia[3];
    float f[4], n3 = 3.141592654f;
    fp = fopen( "data", "w+" );
    /* Write data in:   decimal integer formats
                        decimal, octal, hex integer formats
                        floating-point formats  */
    fprintf( fp, "%d %10d %-10d \n", n1, n1, n1 );
    fprintf( fp, "%d %o %x \n", n2, n2, n2 );
    fprintf( fp, "%f %10.10f %e %5.4e \n", n3, n3, n3, n3 );
    /* Rewind file (to the beginning). */
    rewind( fp );//Moves the file position indicator to the beginning in a file
    /* Read back data. */
    puts( "" );
    fscanf( fp, "%d %d %d", &ia[0], &ia[1], &ia[2] );
    printf( "   %d\t%d\t%d\n", ia[0], ia[1], ia[2] );
    fscanf( fp, "%d %o %x", &ia[0], &ia[1], &ia[2] );
    printf( "   %d\t%d\t%d\n", ia[0], ia[1], ia[2] );
    fscanf( fp, "%f %f %f %f", &f[0], &f[1], &f[2], &f[3] );
    printf( "   %f\t%f\t%f\t%f\n", f[0], f[1], f[2], f[3] );
    fclose( fp );
//----------------------------------------------------------//
    #define SIZE 20
    int n;
    float d[SIZE];
    for( n = 0; n < SIZE; ++n )                 /* Fill array with roots. */
    {
      d[n] = (float)std::sqrt( (double)n );
    }
    for( n = 0; n < SIZE; ++n ) printf( "1:%7.3f ", d[n] );
    fp = fopen( "out1", "wb+" );                 /* Open file. */
    fwrite( d, sizeof( float ), SIZE, fp );     /* Write it to file. */
    //exit(0);
    rewind( fp );                               /* Rewind file. */
    for( n = 0; n < SIZE; ++n ) d[n]=0.0;
    fread( d, sizeof( float ), SIZE, fp );      /* Read back data. */
    for( n = 0; n < SIZE; ++n )                 /* Print array. */
    {
      printf( "2:%d: %7.3f\n", n, d[n] );
    }
    fclose( fp );
    //exit(0);
//----------------------------------------------------------//
    double dd = 12.23;
    int ii = 101;
    long l = 123023L;
    if((fp=fopen("test", "wb+"))==NULL)
    {
      printf("Ошибка при открытии файла.\n");
      exit(1);
    }
    fwrite(&dd, sizeof(double), 1, fp);
    fwrite(&ii, sizeof(int), 1, fp);
    fwrite(&l, sizeof(long), 1, fp);
    rewind(fp);
    dd=0.0;
    ii=0;
    l=0;
    fread(&dd, sizeof(double), 1, fp);
    fread(&ii, sizeof(int), 1, fp);
    fread(&l, sizeof(long), 1, fp);
    printf("%f %d %ld", dd, ii, l);
    fclose(fp);
//----------------------------------------------------------//

    #define MAX 256
    FILE *src, *dst;
    char b[MAX];
    // Try to open source and destination files.
    if ( ( src = fopen( "test", "r" )) == NULL )
    {
      puts( "Can't open input file." );
      exit(1);
    }
    if ( (dst = fopen( "outfile.txt", "w" )) == NULL )
    {
      puts( "Can't open output file." );
      fclose( src );
      exit(1);
    }
    // Copy one file to the next.
    while( ( fgets( b, MAX, src ) ) != NULL )
    {
      fputs( b, dst );
      std::cout<<"fgets()"<<std::endl;
    }
    // All done, close up shop.
    fclose( src );
    fclose( dst );
    //exit(0);
//-----rename file:---------------------------------------------------//
    fp = fopen("data.txt", "w"); // create file "data.txt"
    if(!fp) { perror("data.txt"); return 1; }
    fputc('a', fp); // write a character to "data.txt"
    fclose(fp);
    int rc = rename("data.txt", "to.txt");
    if(rc) { perror("rename"); return 1; }
    fp = fopen("to.txt", "r");
    if(!fp) { perror("to.txt"); return 1; }
    printf("%c\n", fgetc(fp)); // read from "to.txt"
    fclose(fp);
    //exit(0);
//----------------------------------------------------------//

    shop[0].name="apples", shop[1].name="bananas", shop[2].name="bread";
    shop[0].qty=10, shop[1].qty=5, shop[2].qty=1;
    shop[0].price=23.4, shop[1].price=25, shop[2].price=10.3;

    int i=0;
    file = fopen("fprintf.txt", "w");

    std::cout<<"file="<<file<<" "<<scanf("%s%u%f", shop[i].name, &(shop[i].qty), &(shop[i].price))
             <<" "<<shop[i].name<<" "<<" "<<std::endl;

    while (scanf ("%s%u%f", shop[i].name, &(shop[i].qty), &(shop[i].price)) != EOF)
    {
      fprintf(file, "%s %u %.2f\n", shop[i].name, shop[i].qty, shop[i].price);
      i++;
      std::cout<<"i="<<i<<std::endl;
    }

    exit(0);

    file = fopen("fscanf.txt", "r");
    while (fscanf (file, "%s%u%f", shop[i].name, &(shop[i].qty), &(shop[i].price)) != EOF)
    {
      printf("%s %u %.2f\n", shop[i].name, shop[i].qty, shop[i].price);
      i++;
    }
    //
    int N=80;
    char arr[N];
    file = fopen("fscanf.txt", "r");
    while (fgets (arr, N, file) != NULL)
           printf("%s", arr);
    printf("\n");
    fclose(file);
    //
    while ((arr[i] = fgetc (file)) != EOF)
    {
      if (arr[i] == '\n')
      {
        arr[i] = '\0';
         printf("%s\n",arr);
          i = 0;
       }
       else i++;
     }
     arr[i] = '\0';
     printf("%s\n",arr);


    return 0;
}
