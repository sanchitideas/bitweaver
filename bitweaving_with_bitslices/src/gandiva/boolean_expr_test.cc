// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#include <gtest/gtest.h>
#include "arrow/memory_pool.h"
#include "arrow/status.h"
#include <fstream>
#include "gandiva/projector.h"
#include "gandiva/tests/test_util.h"
#include "gandiva/tree_expr_builder.h"

namespace gandiva {

using arrow::boolean;
using arrow::int32;

static inline uint64_t rdtscp( uint32_t & aux ) {
    uint64_t rax,rdx;
    asm volatile ( "rdtscp\n" : "=a" (rax), "=d" (rdx), "=c" (aux) : : );
    return (rdx << 32) + rax;
}

class TestBooleanExpr : public ::testing::Test {
 public:
  void SetUp() { pool_ = arrow::default_memory_pool(); }

 protected:
  arrow::MemoryPool* pool_;
};

TEST_F(TestBooleanExpr, TestSpeed) {
  auto PAGE_SIZE = 4096;
  auto shipdate_field = arrow::field("date", arrow::int16());
  auto discount_field = arrow::field("discount", arrow::int8());
  auto quantity_field = arrow::field("quantity", arrow::int8());
  auto schema = arrow::schema({shipdate_field, discount_field, quantity_field});
  auto lowdate = TreeExprBuilder::MakeLiteral((int16_t)729);
  auto highdate = TreeExprBuilder::MakeLiteral((int16_t)1095);
  auto lowdiscount = TreeExprBuilder::MakeLiteral((int8_t)4);
  auto highdiscount = TreeExprBuilder::MakeLiteral((int8_t)8);
  auto quantity = TreeExprBuilder::MakeLiteral((int8_t)21);


  auto node_shipdate = TreeExprBuilder::MakeField(shipdate_field);
  auto node_discount = TreeExprBuilder::MakeField(discount_field);
  auto node_quantity = TreeExprBuilder::MakeField(quantity_field);

  auto shipdate_gt =
      TreeExprBuilder::MakeFunction("greater_than", {node_shipdate, lowdate}, 
        boolean());
  auto shipdate_lt =
      TreeExprBuilder::MakeFunction("less_than", {node_shipdate, highdate}, 
        boolean());
  auto discount_gt =
      TreeExprBuilder::MakeFunction("greater_than", {node_discount, lowdiscount}, 
        boolean());
  auto discount_lt =
      TreeExprBuilder::MakeFunction("less_than", {node_discount, highdiscount}, 
        boolean());
  auto quantity_lt =
      TreeExprBuilder::MakeFunction("less_than", {node_quantity, quantity}, 
        boolean());

  // output fields
  auto field_result = field("res", boolean());

  auto node_and = TreeExprBuilder::MakeAnd({shipdate_gt, shipdate_lt,
    discount_gt, discount_lt, quantity_lt});
  auto expr = TreeExprBuilder::MakeExpression(node_and, field_result);

  // Build a projector for the expressions.
  std::shared_ptr<Projector> projector;
  auto status = Projector::Make(schema, {expr}, TestConfiguration(), 
    &projector);
  EXPECT_TRUE(status.ok());

  auto num_records = 59986052;
  uint16_t* shipdatePtr = static_cast<uint16_t*>(aligned_alloc(PAGE_SIZE, 
    sizeof(uint16_t) * num_records));
  uint8_t* quantityPtr = static_cast<uint8_t*>(aligned_alloc(PAGE_SIZE, 
    sizeof(uint8_t) * num_records));
  uint8_t* discountPtr = static_cast<uint8_t*>(aligned_alloc(PAGE_SIZE,
    sizeof(uint8_t) * num_records));
  std::vector<int16_t> shipdate_vector(shipdatePtr, shipdatePtr + num_records);
  std::vector<int8_t> discount_vector(discountPtr, discountPtr + num_records);
  std::vector<int8_t> quantity_vector(quantityPtr, quantityPtr + num_records);

  uint8_t* shipdateValidityPtr = static_cast<uint8_t*>(aligned_alloc(PAGE_SIZE, 
    num_records));
  uint8_t* discountValidityPtr = static_cast<uint8_t*>(aligned_alloc(PAGE_SIZE, 
    num_records));
  uint8_t* quantityValidityPtr = static_cast<uint8_t*>(aligned_alloc(PAGE_SIZE, 
    num_records));


  std::ifstream shipdateFile("arrow_shipdate.col", std::ios::binary);
  std::streamsize size = sizeof(uint16_t) * num_records;
  shipdateFile.seekg(0, std::ios::beg);
  shipdateFile.read(static_cast<char*>((void*)shipdatePtr), size);

  std::ifstream discountFile("arrow_discount.col", std::ios::binary);
  size = sizeof(uint8_t) * num_records;
  discountFile.seekg(0, std::ios::beg);
  discountFile.read(static_cast<char*>((void*)discountPtr), size);

  std::ifstream quantityFile("arrow_quantity.col", std::ios::binary);
  size = sizeof(uint8_t) * num_records;
  quantityFile.seekg(0, std::ios::beg);
  quantityFile.read(static_cast<char*>((void*)quantityPtr), size);


  std::vector<bool> shipdate_validity_vector(shipdateValidityPtr, 
    shipdateValidityPtr + num_records);
  std::vector<bool> discount_validity_vector(discountValidityPtr, 
    discountValidityPtr + num_records);
  std::vector<bool> quantity_validity_vector(quantityValidityPtr, 
    quantityValidityPtr + num_records);
  std::fill (shipdate_validity_vector.begin(), shipdate_validity_vector.end(), 
    true);
  std::fill (discount_validity_vector.begin(), discount_validity_vector.end(), 
    true);
  std::fill (quantity_validity_vector.begin(), quantity_validity_vector.end(), 
    true);

  auto shipdate_array = MakeArrowArrayInt16(shipdate_vector, 
    shipdate_validity_vector);
  auto discount_array = MakeArrowArrayInt8(discount_vector, 
    discount_validity_vector);
  auto quantity_array = MakeArrowArrayInt8(quantity_vector, 
    quantity_validity_vector);
  auto in_batch = arrow::RecordBatch::Make(schema, num_records, 
    {shipdate_array, discount_array, quantity_array});
 
  arrow::ArrayVector outputs;
  uint32_t placeHolder = 0;
  auto startTime = std::chrono::high_resolution_clock::now();
  uint64_t start = rdtscp(placeHolder);
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  uint64_t end = rdtscp(placeHolder);
  auto stopTime = std::chrono::high_resolution_clock::now();
  uint64_t diff = (uint64_t)(end - start);
  std::cout << "\n\nIn AVX-512 based Arrow"
    " The number of cycles are " << diff << " and the number of cycles "
    "per code are " << diff/(num_records * 1.0) << std::endl;
  auto duration = 
    std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    std::cout << "\n\nTime taken by AVX-512 Arrow is " <<
    duration.count() << " microseconds\n\n" << std::endl;
  EXPECT_TRUE(status.ok());
}




TEST_F(TestBooleanExpr, SimpleAnd) {
  // schema for input fields
  auto fielda = field("a", int32());
  auto fieldb = field("b", int32());
  auto schema = arrow::schema({fielda, fieldb});




  // output fields
  auto field_result = field("res", boolean());

  // build expression.
  // (a > 0) && (b > 0)
  auto node_a = TreeExprBuilder::MakeField(fielda);
  auto node_b = TreeExprBuilder::MakeField(fieldb);
  auto literal_0 = TreeExprBuilder::MakeLiteral((int32_t)0);
  auto a_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, literal_0}, boolean());
  auto b_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_b, literal_0}, boolean());

  auto node_and = TreeExprBuilder::MakeAnd({a_gt_0, b_gt_0});
  auto expr = TreeExprBuilder::MakeExpression(node_and, field_result);

  // Build a projector for the expressions.
  std::shared_ptr<Projector> projector;
  auto status = Projector::Make(schema, {expr}, TestConfiguration(), &projector);
  EXPECT_TRUE(status.ok());

  // FALSE_VALID && ?  => FALSE_VALID
  int num_records = 4;
  auto arraya = MakeArrowArrayInt32({-2, -2, -2, -2}, {true, true, true, true});
  auto arrayb = MakeArrowArrayInt32({-2, -2, 2, 2}, {true, false, true, false});
  auto exp = MakeArrowArrayBool({false, false, false, false}, {true, true, true, true});
  auto in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});

  arrow::ArrayVector outputs;
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));

  // FALSE_INVALID && ?
  num_records = 4;
  arraya = MakeArrowArrayInt32({-2, -2, -2, -2}, {false, false, false, false});
  arrayb = MakeArrowArrayInt32({-2, -2, 2, 2}, {true, false, true, false});
  exp = MakeArrowArrayBool({false, false, false, false}, {true, false, false, false});
  in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});
  outputs.clear();
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));

  // TRUE_VALID && ?
  num_records = 4;
  arraya = MakeArrowArrayInt32({2, 2, 2, 2}, {true, true, true, true});
  arrayb = MakeArrowArrayInt32({-2, -2, 2, 2}, {true, false, true, false});
  exp = MakeArrowArrayBool({false, false, true, false}, {true, false, true, false});
  in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});
  outputs.clear();
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));

  // TRUE_INVALID && ?
  num_records = 4;
  arraya = MakeArrowArrayInt32({2, 2, 2, 2}, {false, false, false, false});
  arrayb = MakeArrowArrayInt32({-2, -2, 2, 2}, {true, false, true, false});
  exp = MakeArrowArrayBool({false, false, false, false}, {true, false, false, false});
  in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});
  outputs.clear();
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));
}

TEST_F(TestBooleanExpr, SimpleOr) {
  // schema for input fields
  auto fielda = field("a", int32());
  auto fieldb = field("b", int32());
  auto schema = arrow::schema({fielda, fieldb});

  // output fields
  auto field_result = field("res", boolean());

  // build expression.
  // (a > 0) || (b > 0)
  auto node_a = TreeExprBuilder::MakeField(fielda);
  auto node_b = TreeExprBuilder::MakeField(fieldb);
  auto literal_0 = TreeExprBuilder::MakeLiteral((int32_t)0);
  auto a_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, literal_0}, boolean());
  auto b_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_b, literal_0}, boolean());

  auto node_or = TreeExprBuilder::MakeOr({a_gt_0, b_gt_0});
  auto expr = TreeExprBuilder::MakeExpression(node_or, field_result);

  // Build a projector for the expressions.
  std::shared_ptr<Projector> projector;
  auto status = Projector::Make(schema, {expr}, TestConfiguration(), &projector);
  EXPECT_TRUE(status.ok());

  // TRUE_VALID && ?  => TRUE_VALID
  int num_records = 4;
  auto arraya = MakeArrowArrayInt32({2, 2, 2, 2}, {true, true, true, true});
  auto arrayb = MakeArrowArrayInt32({-2, -2, 2, 2}, {true, false, true, false});
  auto exp = MakeArrowArrayBool({true, true, true, true}, {true, true, true, true});
  auto in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});

  arrow::ArrayVector outputs;
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));

  // TRUE_INVALID && ?
  num_records = 4;
  arraya = MakeArrowArrayInt32({2, 2, 2, 2}, {false, false, false, false});
  arrayb = MakeArrowArrayInt32({-2, -2, 2, 2}, {true, false, true, false});
  exp = MakeArrowArrayBool({false, false, true, false}, {false, false, true, false});
  in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});
  outputs.clear();
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));

  // FALSE_VALID && ?
  num_records = 4;
  arraya = MakeArrowArrayInt32({-2, -2, -2, -2}, {true, true, true, true});
  arrayb = MakeArrowArrayInt32({-2, -2, 2, 2}, {true, false, true, false});
  exp = MakeArrowArrayBool({false, false, true, false}, {true, false, true, false});
  in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});
  outputs.clear();
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));

  // FALSE_INVALID && ?
  num_records = 4;
  arraya = MakeArrowArrayInt32({-2, -2, -2, -2}, {false, false, false, false});
  arrayb = MakeArrowArrayInt32({-2, -2, 2, 2}, {true, false, true, false});
  exp = MakeArrowArrayBool({false, false, true, false}, {false, false, true, false});
  in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});
  outputs.clear();
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));
}

TEST_F(TestBooleanExpr, AndThree) {
  // schema for input fields
  auto fielda = field("a", int32());
  auto fieldb = field("b", int32());
  auto fieldc = field("c", int32());
  auto schema = arrow::schema({fielda, fieldb, fieldc});

  // output fields
  auto field_result = field("res", boolean());

  // build expression.
  // (a > 0) && (b > 0) && (c > 0)
  auto node_a = TreeExprBuilder::MakeField(fielda);
  auto node_b = TreeExprBuilder::MakeField(fieldb);
  auto node_c = TreeExprBuilder::MakeField(fieldc);
  auto literal_0 = TreeExprBuilder::MakeLiteral((int32_t)0);
  auto a_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, literal_0}, boolean());
  auto b_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_b, literal_0}, boolean());
  auto c_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_c, literal_0}, boolean());

  auto node_and = TreeExprBuilder::MakeAnd({a_gt_0, b_gt_0, c_gt_0});
  auto expr = TreeExprBuilder::MakeExpression(node_and, field_result);

  // Build a projector for the expressions.
  std::shared_ptr<Projector> projector;
  auto status = Projector::Make(schema, {expr}, TestConfiguration(), &projector);
  EXPECT_TRUE(status.ok());

  int num_records = 8;
  std::vector<bool> validity({true, true, true, true, true, true, true, true});
  auto arraya = MakeArrowArrayInt32({2, 2, 2, 0, 2, 0, 0, 0}, validity);
  auto arrayb = MakeArrowArrayInt32({2, 2, 0, 2, 0, 2, 0, 0}, validity);
  auto arrayc = MakeArrowArrayInt32({2, 0, 2, 2, 0, 0, 2, 0}, validity);
  auto exp = MakeArrowArrayBool({true, false, false, false, false, false, false, false},
                                validity);

  auto in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb, arrayc});

  arrow::ArrayVector outputs;
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));
}

TEST_F(TestBooleanExpr, OrThree) {
  // schema for input fields
  auto fielda = field("a", int32());
  auto fieldb = field("b", int32());
  auto fieldc = field("c", int32());
  auto schema = arrow::schema({fielda, fieldb, fieldc});

  // output fields
  auto field_result = field("res", boolean());

  // build expression.
  // (a > 0) || (b > 0) || (c > 0)
  auto node_a = TreeExprBuilder::MakeField(fielda);
  auto node_b = TreeExprBuilder::MakeField(fieldb);
  auto node_c = TreeExprBuilder::MakeField(fieldc);
  auto literal_0 = TreeExprBuilder::MakeLiteral((int32_t)0);
  auto a_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, literal_0}, boolean());
  auto b_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_b, literal_0}, boolean());
  auto c_gt_0 =
      TreeExprBuilder::MakeFunction("greater_than", {node_c, literal_0}, boolean());

  auto node_or = TreeExprBuilder::MakeOr({a_gt_0, b_gt_0, c_gt_0});
  auto expr = TreeExprBuilder::MakeExpression(node_or, field_result);

  // Build a projector for the expressions.
  std::shared_ptr<Projector> projector;
  auto status = Projector::Make(schema, {expr}, TestConfiguration(), &projector);
  EXPECT_TRUE(status.ok());

  int num_records = 8;
  std::vector<bool> validity({true, true, true, true, true, true, true, true});
  auto arraya = MakeArrowArrayInt32({2, 2, 2, 0, 2, 0, 0, 0}, validity);
  auto arrayb = MakeArrowArrayInt32({2, 2, 0, 2, 0, 2, 0, 0}, validity);
  auto arrayc = MakeArrowArrayInt32({2, 0, 2, 2, 0, 0, 2, 0}, validity);
  auto exp =
      MakeArrowArrayBool({true, true, true, true, true, true, true, false}, validity);

  auto in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb, arrayc});

  arrow::ArrayVector outputs;
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));
}

TEST_F(TestBooleanExpr, BooleanAndInsideIf) {
  // schema for input fields
  auto fielda = field("a", int32());
  auto fieldb = field("b", int32());
  auto schema = arrow::schema({fielda, fieldb});

  // output fields
  auto field_result = field("res", boolean());

  // build expression.
  // if (a > 2 && b > 2)
  //   a > 3 && b > 3
  // else
  //   a > 1 && b > 1
  auto node_a = TreeExprBuilder::MakeField(fielda);
  auto node_b = TreeExprBuilder::MakeField(fieldb);
  auto literal_1 = TreeExprBuilder::MakeLiteral((int32_t)1);
  auto literal_2 = TreeExprBuilder::MakeLiteral((int32_t)2);
  auto literal_3 = TreeExprBuilder::MakeLiteral((int32_t)3);
  auto a_gt_1 =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, literal_1}, boolean());
  auto a_gt_2 =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, literal_2}, boolean());
  auto a_gt_3 =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, literal_3}, boolean());
  auto b_gt_1 =
      TreeExprBuilder::MakeFunction("greater_than", {node_b, literal_1}, boolean());
  auto b_gt_2 =
      TreeExprBuilder::MakeFunction("greater_than", {node_b, literal_2}, boolean());
  auto b_gt_3 =
      TreeExprBuilder::MakeFunction("greater_than", {node_b, literal_3}, boolean());

  auto and_1 = TreeExprBuilder::MakeAnd({a_gt_1, b_gt_1});
  auto and_2 = TreeExprBuilder::MakeAnd({a_gt_2, b_gt_2});
  auto and_3 = TreeExprBuilder::MakeAnd({a_gt_3, b_gt_3});

  auto node_if = TreeExprBuilder::MakeIf(and_2, and_3, and_1, arrow::boolean());
  auto expr = TreeExprBuilder::MakeExpression(node_if, field_result);

  // Build a projector for the expressions.
  std::shared_ptr<Projector> projector;
  auto status = Projector::Make(schema, {expr}, TestConfiguration(), &projector);
  EXPECT_TRUE(status.ok());

  int num_records = 4;
  std::vector<bool> validity({true, true, true, true});
  auto arraya = MakeArrowArrayInt32({4, 4, 2, 1}, validity);
  auto arrayb = MakeArrowArrayInt32({5, 3, 3, 1}, validity);
  auto exp = MakeArrowArrayBool({true, false, true, false}, validity);

  auto in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});

  arrow::ArrayVector outputs;
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));
}

TEST_F(TestBooleanExpr, IfInsideBooleanAnd) {
  // schema for input fields
  auto fielda = field("a", int32());
  auto fieldb = field("b", int32());
  auto schema = arrow::schema({fielda, fieldb});

  // output fields
  auto field_result = field("res", boolean());

  // build expression.
  // (if (a > b) a > 3 else b > 3) && (if (a > b) a > 2 else b > 2)

  auto node_a = TreeExprBuilder::MakeField(fielda);
  auto node_b = TreeExprBuilder::MakeField(fieldb);
  auto literal_2 = TreeExprBuilder::MakeLiteral((int32_t)2);
  auto literal_3 = TreeExprBuilder::MakeLiteral((int32_t)3);
  auto a_gt_b =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, node_b}, boolean());
  auto a_gt_2 =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, literal_2}, boolean());
  auto a_gt_3 =
      TreeExprBuilder::MakeFunction("greater_than", {node_a, literal_3}, boolean());
  auto b_gt_2 =
      TreeExprBuilder::MakeFunction("greater_than", {node_b, literal_2}, boolean());
  auto b_gt_3 =
      TreeExprBuilder::MakeFunction("greater_than", {node_b, literal_3}, boolean());

  auto if_3 = TreeExprBuilder::MakeIf(a_gt_b, a_gt_3, b_gt_3, arrow::boolean());
  auto if_2 = TreeExprBuilder::MakeIf(a_gt_b, a_gt_2, b_gt_2, arrow::boolean());
  auto node_and = TreeExprBuilder::MakeAnd({if_3, if_2});
  auto expr = TreeExprBuilder::MakeExpression(node_and, field_result);

  // Build a projector for the expressions.
  std::shared_ptr<Projector> projector;
  auto status = Projector::Make(schema, {expr}, TestConfiguration(), &projector);
  EXPECT_TRUE(status.ok());

  int num_records = 4;
  std::vector<bool> validity({true, true, true, true});
  auto arraya = MakeArrowArrayInt32({4, 3, 3, 2}, validity);
  auto arrayb = MakeArrowArrayInt32({3, 4, 2, 3}, validity);
  auto exp = MakeArrowArrayBool({true, true, false, false}, validity);

  auto in_batch = arrow::RecordBatch::Make(schema, num_records, {arraya, arrayb});

  arrow::ArrayVector outputs;
  status = projector->Evaluate(*in_batch, pool_, &outputs);
  EXPECT_TRUE(status.ok());
  EXPECT_ARROW_ARRAY_EQUALS(exp, outputs.at(0));
}

}  // namespace gandiva
